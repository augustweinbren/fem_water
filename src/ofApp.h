#pragma once

#include "ofMain.h"

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
    float waveSpeed = 1.0f;
    float damping = 0.999f;
    
    // Adjust this to start with a larger wave:
    float initialWaveAmplitude = 20.0f;
    
    // Adjust this to control random noise impulses:
    float noiseProbability     = 0.005f;
    float noiseMagnitude       = 20.0f;
    
    // Finite Element data
    vector<ofVec2f> nodes;       // Node positions
    vector<float>   u;           // Displacement at each node
    vector<float>   uOld;        // Displacement at previous time step
    vector<float>   uVel;        // Velocity at each node
    vector<Triangle> triangles;  // Triangular elements
    
    // Utility methods
    void buildMesh(int gridWidth, int gridHeight, float spacing);
    void applyInitialConditions();
    void solveWaveEquationFEM();
    inline float cross2D(const ofVec2f &v1, const ofVec2f &v2){
        return v1.x * v2.y - v1.y * v2.x;
    }
};
