#pragma once

#include "ofMain.h"

// Include Eigen (make sure you have the Eigen library installed and accessible)
#include <Eigen/Dense>
#include <Eigen/Sparse>

struct Triangle {
    int a, b, c;
};

class ofApp : public ofBaseApp{

	public:
		void setup() override;
		void update() override;
		void draw() override;
		void exit() override;

		void keyPressed(int key) override;
		void keyReleased(int key) override;
		void mouseMoved(int x, int y ) override;
		void mouseDragged(int x, int y, int button) override;
		void mousePressed(int x, int y, int button) override;
		void mouseReleased(int x, int y, int button) override;
		void mouseScrolled(int x, int y, float scrollX, float scrollY) override;
		void mouseEntered(int x, int y) override;
		void mouseExited(int x, int y) override;
		void windowResized(int w, int h) override;
		void dragEvent(ofDragInfo dragInfo) override;
		void gotMessage(ofMessage msg) override;
		
private:
    float timeStep = 0.01f;
    float waveSpeed = 10.0f;
    float damping = 0.9999f;
    
    // Adjust this to start with a larger wave:
    float initialWaveAmplitude = 100.0f;
    
    // Adjust this to control random noise impulses:
    float noiseProbability     = 0.01f;
    float noiseMagnitude       = 20.0f;
    
    // Finite Element data
    vector<ofVec2f> nodes;       // Node positions
    vector<Triangle> triangles;  // Triangular elements
    
    // Displacement (U), Velocity (V), Acceleration (A)
    // stored as 1D arrays (one scalar DOF per node) for wave simulation
    vector<float> U;
    vector<float> V;
    vector<float> A;
    
    // Lumped mass (M[i]) for each node i
    vector<float> M;
    
    // Sparse stiffness matrix (K)
    // We store only the upper-tri or full. Here we use a simple triplet-based approach for assembly.
    Eigen::SparseMatrix<float> K;
    
    // Utility methods
    void buildMesh(int gridWidth, int gridHeight, float spacing);
    
    // Assemble M (lumped) and K (sparse) for the 2D wave equation
    void assembleSystem();
    
    //    No longer needed (?)
    void applyInitialConditions();
    
    void solveWaveEquationFEM();
    
    inline float cross2D(const ofVec2f &v1, const ofVec2f &v2){
        return v1.x * v2.y - v1.y * v2.x;
    }
};
