#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
    // Set up an orthographic camera for 2D
    ofSetVerticalSync(true);
    ofBackground(0);

    // Build a small mesh of, say, 20x20 nodes for demonstration
    buildMesh(20, 20, 10.0f);

    // Apply initial displacement/velocity
    applyInitialConditions();
}

//--------------------------------------------------------------
void ofApp::update(){
    // Update the wave equation using a simplified FEM approach
    solveWaveEquationFEM();
    
    // Add random noise impulses from outside
    if(ofRandom(1.0f) < noiseProbability) {
        int idx = (int)ofRandom(nodes.size());
        u[idx] += ofRandom(-noiseMagnitude, noiseMagnitude);
    }
}

//--------------------------------------------------------------
void ofApp::draw(){
    // Translate to center
    ofPushMatrix();
    ofTranslate(ofGetWidth()/2, ofGetHeight()/2);

    // Draw the triangular mesh
    ofSetColor(0, 150, 255);
    ofNoFill();
    for(const auto & tri : triangles) {
        ofDrawTriangle(nodes[tri.a].x, nodes[tri.a].y,
                       nodes[tri.b].x, nodes[tri.b].y,
                       nodes[tri.c].x, nodes[tri.c].y);
    }

    // Draw the displacement as circles at each node
    ofFill();
    for(int i=0; i<nodes.size(); i++){
        // Map displacement to a brightness value or vertical offset
        float brightness = ofMap(u[i], -2.0f, 2.0f, 50, 255, true);
        ofSetColor(brightness, brightness, 255);
        ofDrawCircle(nodes[i].x, nodes[i].y - u[i], 2);
    }

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

    // Create triangle connectivity (two triangles per cell)
    for(int y=0; y<gridHeight-1; y++){
        for(int x=0; x<gridWidth-1; x++){
            int i0 = y     * gridWidth + x;
            int i1 = y     * gridWidth + (x+1);
            int i2 = (y+1) * gridWidth + x;
            int i3 = (y+1) * gridWidth + (x+1);

            // Triangle 1
            Triangle t1;
            t1.a = i0; t1.b = i1; t1.c = i2;
            triangles.push_back(t1);

            // Triangle 2
            Triangle t2;
            t2.a = i1; t2.b = i3; t2.c = i2;
            triangles.push_back(t2);
        }
    }

    // Initialize displacement arrays
    u.assign(nodes.size(), 0.0f);
    uOld.assign(nodes.size(), 0.0f);
    uVel.assign(nodes.size(), 0.0f);
}

//--------------------------------------------------------------
void ofApp::applyInitialConditions(){
    // Simple bump in the center
    if(!nodes.empty()) {
        int center = nodes.size() / 2;
        u[center] = initialWaveAmplitude;  // a small initial displacement
    }
}

//--------------------------------------------------------------
// Extremely simplified explicit FEM approach for a 2D wave
// (for demonstration only, not a fully correct or stable solver)
void ofApp::solveWaveEquationFEM(){
    // We do a naive explicit approach:
    //    uAcc ~ waveSpeed^2 * Laplacian(u)
    // with some damping and a simple free boundary condition.

    vector<float> uAcc(nodes.size(), 0.0f);

    // For each triangle, approximate Laplacian contribution
    // using a simple linear approach (not full detail).
    for(const auto & tri : triangles){
        // Indices
        int iA = tri.a;
        int iB = tri.b;
        int iC = tri.c;

        // Positions
        ofVec2f &pA = nodes[iA];
        ofVec2f &pB = nodes[iB];
        ofVec2f &pC = nodes[iC];

        // Displacements
        float uA = u[iA];
        float uB = u[iB];
        float uC = u[iC];

        // Triangle area (2D cross)
        float area = fabs(cross2D(pB - pA, pC - pA)) * 0.5f;

        // Simple finite difference of "neighbors" to approximate
        // partial Laplacian. This is *very* simplified:
        float avgU = (uA + uB + uC) / 3.0f;
        float dA = uA - avgU;
        float dB = uB - avgU;
        float dC = uC - avgU;

        // Contribute to each node's acceleration
        // Weighted by area as a placeholder for actual shape function integration
        uAcc[iA] += -dA * area;
        uAcc[iB] += -dB * area;
        uAcc[iC] += -dC * area;
    }

    // Integrate in time (explicit)
    for(int i=0; i<nodes.size(); i++){
        // Acceleration from approximate Laplacian
        float acc = waveSpeed * waveSpeed * uAcc[i];

        // Simple velocity Verlet or forward Euler update
        uVel[i] = damping * (uVel[i] + acc * timeStep);
        u[i]    = u[i] + uVel[i] * timeStep;
    }

    // Optionally fix outer boundary by clamping displacement
    // Example: clamp edges to zero
    // (In a real FEM solution, we'd handle boundary conditions more rigorously)
    // For demonstration, let it be free or selectively clamp corners:
    // uVel[...] = 0.0f, u[...] = 0.0f if desired.
}

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
