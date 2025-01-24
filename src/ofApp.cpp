#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
    // Set up an orthographic camera for 2D
    ofSetVerticalSync(true);
    ofBackground(0);

    // Build a small mesh of, say, 20x20 nodes for demonstration
    buildMesh(70, 70, 10.0f);

    // Assemble the global mass (lumped) and stiffness matrices
    assembleSystem();

    // Initialize displacements and velocities
    U.resize(nodes.size(), 0.0f);
    V.resize(nodes.size(), 0.0f);
    A.resize(nodes.size(), 0.0f);

    // Give an initial "bump" in the middle
    int center = (int)nodes.size() / 2;
    if(!nodes.empty()) {
        U[center] = initialWaveAmplitude;
    }
}

//--------------------------------------------------------------
void ofApp::update(){
    // Update the wave equation using a simplified FEM approach
    solveWaveEquationFEM();
    
    // Add random noise impulses from outside
    if(ofRandom(1.0f) < noiseProbability) {
        int idx = (int)ofRandom(nodes.size());
        U[idx] *= ofRandom(-noiseMagnitude, noiseMagnitude);
    }
}

//--------------------------------------------------------------
void ofApp::draw(){
    // Translate to center
    ofPushMatrix();
    ofTranslate(ofGetWidth() * 0.5f, ofGetHeight() * 0.5f);

    // Draw the triangular mesh
    ofSetColor(0, 150, 255);
    ofNoFill();
    for(const auto & tri : triangles) {
        auto &pA = nodes[tri.a];
        auto &pB = nodes[tri.b];
        auto &pC = nodes[tri.c];
        ofDrawTriangle(pA.x, pA.y - U[tri.a],
                       pB.x, pB.y - U[tri.b],
                       pC.x, pC.y - U[tri.c]);
    }
//NOTE: Maybe readd the below
//    // Draw the displacement as circles at each node
//    ofFill();
//    for(int i=0; i<nodes.size(); i++){
//        // Map displacement to a brightness value or vertical offset
//        float brightness = ofMap(U[i], -2.0f, 2.0f, 50, 255, true);
//        ofSetColor(brightness, brightness, 255);
//        ofDrawCircle(nodes[i].x, nodes[i].y - U[i], 2);
//    }

    ofPopMatrix();
}

//--------------------------------------------------------------
void ofApp::buildMesh(int gridWidth, int gridHeight, float spacing){
    nodes.clear();
    triangles.clear();

    // Create grid nodes
    for(int y=0; y<gridHeight; y++){
        for(int x=0; x<gridWidth; x++){
            float px = (x - gridWidth * 0.5f) * spacing;
            float py = (y - gridHeight * 0.5f) * spacing;
            nodes.push_back(ofVec2f(px, py));
        }
    }

    // Create triangular elements (two per cell)
    for(int y=0; y<gridHeight-1; y++){
        for(int x=0; x<gridWidth-1; x++){
            int i0 =  y      * gridWidth + x;
            int i1 =  y      * gridWidth + (x+1);
            int i2 = (y+1)   * gridWidth + x;
            int i3 = (y+1)   * gridWidth + (x+1);

            triangles.push_back({i0, i1, i2});
            triangles.push_back({i1, i3, i2});
        }
    }
}

//--------------------------------------------------------------
// Assemble the global stiffness (K) and lumped mass (M)
void ofApp::assembleSystem(){
    // Resize M for # of nodes and initialize to 0
    M.assign(nodes.size(), 0.0f);

    // We'll assemble K using triplets, then build a sparse matrix
    vector<Eigen::Triplet<float>> triplets;
    triplets.reserve(triangles.size() * 9); // each triangle contributes up to 9 entries

    // For each triangle, compute local mass and stiffness
    for(auto &tri : triangles){
        int iA = tri.a;
        int iB = tri.b;
        int iC = tri.c;

        // Node coords
        auto &pA = nodes[iA];
        auto &pB = nodes[iB];
        auto &pC = nodes[iC];

        // 2*Area (signed)
        float twiceAreaSigned = cross2D((pB - pA), (pC - pA));
        float area = 0.5f * fabs(twiceAreaSigned);

        // Compute b_i, c_i for the standard linear shape gradients
        // b1 = yB - yC, b2 = yC - yA, b3 = yA - yB
        float b1 = pB.y - pC.y;
        float b2 = pC.y - pA.y;
        float b3 = pA.y - pB.y;

        // c1 = xC - xB, c2 = xA - xC, c3 = xB - xA
        float c1 = pC.x - pB.x;
        float c2 = pA.x - pC.x;
        float c3 = pB.x - pA.x;

        // Local stiffness matrix K_e for the Laplacian:
        // K_e = (1 / (2*Area)) * [ [b1b1 + c1c1, b1b2 + c1c2, ...], ... ]
        // We then multiply by waveSpeed^2 after assembly if desired,
        // but let's incorporate waveSpeed^2 directly here:
        float factor = waveSpeed * waveSpeed / (2.0f * area);

        // Each entry K_ij = factor * (b_i*b_j + c_i*c_j)
        float k11 = factor * (b1*b1 + c1*c1);
        float k12 = factor * (b1*b2 + c1*c2);
        float k13 = factor * (b1*b3 + c1*c3);
        float k22 = factor * (b2*b2 + c2*c2);
        float k23 = factor * (b2*b3 + c2*c3);
        float k33 = factor * (b3*b3 + c3*c3);

        // Add to triplet list (symmetric matrix)
        auto addSym = [&](int r, int c, float val){
            triplets.push_back(Eigen::Triplet<float>(r, c, val));
            if(r != c){
                triplets.push_back(Eigen::Triplet<float>(c, r, val));
            }
        };

        addSym(iA, iA, k11);
        addSym(iA, iB, k12);
        addSym(iA, iC, k13);
        addSym(iB, iB, k22);
        addSym(iB, iC, k23);
        addSym(iC, iC, k33);

        // Local mass matrix M_e (consistent):
        // M_e = (rho * thickness * area) / 12 * [ [2,1,1],[1,2,1],[1,1,2] ]
        // For simplicity, assume density=1, thickness=1 => M_e = area/12 * ...
        // We'll do a LUMPED mass approach, so each node gets 1/3 of the total element mass.
        // Total element mass = area * (1). So each of the 3 nodes gets (area/3).
        float lumpedVal = area / 3.0f;
        M[iA] += lumpedVal;
        M[iB] += lumpedVal;
        M[iC] += lumpedVal;
    }

    // Build K from triplets
    K.resize(nodes.size(), nodes.size());
    K.setFromTriplets(triplets.begin(), triplets.end());
    K.makeCompressed();
}

//--------------------------------------------------------------
// Explicit time integration with lumped mass
void ofApp::solveWaveEquationFEM(){
    // First compute acceleration A = -M^-1 K U
    // Since M is lumped (diagonal), M^-1 is simply 1/M[i].
    // We'll store the result in A.
    // Then V += A * dt, and U += V * dt, plus damping.
    // This is a standard explicit integrator for the wave equation.

    // Multiply K * U into a temp vector (K is NxN, U is Nx1)
    Eigen::VectorXf vecU(nodes.size()), vecKU(nodes.size());
    for(size_t i=0; i<nodes.size(); i++){
        vecU(i) = U[i];
    }

    vecKU = K * vecU;  // matrix-vector product

    // Compute acceleration A[i] = -(1 / M[i]) * vecKU[i]
    for(size_t i=0; i<nodes.size(); i++){
        A[i] = - (1.0f / M[i]) * vecKU((int)i);
    }

    // Update velocity and displacement
    for(size_t i=0; i<nodes.size(); i++){
        V[i] = damping * (V[i] + A[i] * timeStep);
        U[i] = U[i] + V[i] * timeStep;
    }
}

//--------------------------------------------------------------
//void ofApp::applyInitialConditions(){
//    // Simple bump in the center
//    if(!nodes.empty()) {
//        int center = nodes.size() / 2;
//        u[center] = initialWaveAmplitude;  // a small initial displacement
//    }
//}

//--------------------------------------------------------------
// Extremely simplified explicit FEM approach for a 2D wave
// (for demonstration only, not a fully correct or stable solver)
//void ofApp::solveWaveEquationFEM(){
//    // We do a naive explicit approach:
//    //    uAcc ~ waveSpeed^2 * Laplacian(u)
//    // with some damping and a simple free boundary condition.
//
//    vector<float> uAcc(nodes.size(), 0.0f);
//
//    // For each triangle, approximate Laplacian contribution
//    // using a simple linear approach (not full detail).
//    for(const auto & tri : triangles){
//        // Indices
//        int iA = tri.a;
//        int iB = tri.b;
//        int iC = tri.c;
//
//        // Positions
//        ofVec2f &pA = nodes[iA];
//        ofVec2f &pB = nodes[iB];
//        ofVec2f &pC = nodes[iC];
//
//        // Displacements
//        float uA = u[iA];
//        float uB = u[iB];
//        float uC = u[iC];
//
//        // Triangle area (2D cross)
//        float area = fabs(cross2D(pB - pA, pC - pA)) * 0.5f;
//
//        // Simple finite difference of "neighbors" to approximate
//        // partial Laplacian. This is *very* simplified:
//        float avgU = (uA + uB + uC) / 3.0f;
//        float dA = uA - avgU;
//        float dB = uB - avgU;
//        float dC = uC - avgU;
//
//        // Contribute to each node's acceleration
//        // Weighted by area as a placeholder for actual shape function integration
//        uAcc[iA] += -dA * area;
//        uAcc[iB] += -dB * area;
//        uAcc[iC] += -dC * area;
//    }
//
//    // Integrate in time (explicit)
//    for(int i=0; i<nodes.size(); i++){
//        // Acceleration from approximate Laplacian
//        float acc = waveSpeed * waveSpeed * uAcc[i];
//
//        // Simple velocity Verlet or forward Euler update
//        uVel[i] = damping * (uVel[i] + acc * timeStep);
//        u[i]    = u[i] + uVel[i] * timeStep;
//    }
//
//    // Optionally fix outer boundary by clamping displacement
//    // Example: clamp edges to zero
//    // (In a real FEM solution, we'd handle boundary conditions more rigorously)
//    // For demonstration, let it be free or selectively clamp corners:
//    // uVel[...] = 0.0f, u[...] = 0.0f if desired.
//}

//--------------------------------------------------------------
void ofApp::exit(){

}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseScrolled(int x, int y, float scrollX, float scrollY){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
