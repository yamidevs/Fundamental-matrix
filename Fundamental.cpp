// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse
// Date:     2013/10/08

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}


Matrix<float> compute8pointAlgorithm(vector<Match>& matches){

    size_t size = matches.size() < 9 ? 9 : matches.size();
    Matrix<float> A(size, 9);

    float norm = 0.0001f;

    Matrix<float> N(3,3);
    N = N.Zero(3,3);
    N(0,0) = norm;
    N(1,1) = norm;
    N(2,2) = 1;

    int i = 0;
    // compute matrix A = x x'
    for(auto it = matches.begin();it != matches.end();it++,i++){

         float x,xp,y,yp;
         // normalize our points to avoid numerical problem
         x  = it->x1 * norm ;
         xp = it->x2 * norm;
         y  = it->y1 * norm;
         yp = it->y2 * norm;

         A(i,0) = x * xp;
         A(i,1) = x * yp;
         A(i,2) = x;
         A(i,3) = y * xp;
         A(i,4) = y * yp;
         A(i,5) = y;
         A(i,6) = xp;
         A(i,7) = yp;
         A(i,8) = 1;
    }

    if(matches.size()<9){
        // add 9 equation to be square matrix 9x9
        for(int i=0; i<9; i++){
            A(8,i) = 0;
        }
    }

    Matrix<float> U,V;
    Vector<float> Sigma;

    svd(A,U,Sigma,V);

    Matrix<float> F(3,3);
    // extract F from V ( get our matrix F from the last row V)
    F(0,0) = V.getRow(V.nrow() - 1)[0];
    F(0,1) = V.getRow(V.nrow() - 1)[1];
    F(0,2) = V.getRow(V.nrow() - 1)[2];
    F(1,0) = V.getRow(V.nrow() - 1)[3];
    F(1,1) = V.getRow(V.nrow() - 1)[4];
    F(1,2) = V.getRow(V.nrow() - 1)[5];
    F(2,0) = V.getRow(V.nrow() - 1)[6];
    F(2,1) = V.getRow(V.nrow() - 1)[7];
    F(2,2) = V.getRow(V.nrow() - 1)[8];

    svd(F,U,Sigma,V);
    // force F to be rank 2
    Sigma[2] = 0;
    // recompute F
    F = U * Diagonal(Sigma) * V;

    //denormalize
    F = transpose(N) * F * N;
    return F;

}
// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
Matrix<float> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    Matrix<float> bestF;
    vector<Match> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS
    cout <<"test : "<< matches.size() << endl;
    int n = 0;
    while(n < Niter){

        vector<Match> inliers;
        // take 8 element from matches to compute 8point Algorithm
        for(int i = 0;i<8;i++){
            int randomIndex = rand() % matches.size();
            inliers.push_back(matches[randomIndex]);
        }
        //compute matrix F
        Matrix<float>  F = compute8pointAlgorithm(inliers);
        inliers.clear();
        for(auto it = matches.begin();it != matches.end();it++){
             Vector<float> x(3);
             x[0] = it->x1;
             x[1] = it->y1;
             x[2] = 1;

             Vector<float> xp(3);
             xp[0] = it->x2;
             xp[1] = it->y2;
             xp[2] = 1;

             xp = F * xp;
             // compute the distance between x and Fx' and should be less than distMax , so if it's true,it should be an inliers
             if((pow((x*xp),2)/(pow(xp[0],2)+pow(xp[1],2))) <= pow(distMax,2)){
                 inliers.push_back(*it);
             }
        }
        // each time we found a better model we recompute
        if(inliers.size() > bestInliers.size()){
             bestInliers = inliers;
             float ratio = inliers.size()/(float)matches.size();
             Niter = (int)log(BETA)/log(1-pow(ratio, 8));
        }

        n++;
    }

    // recompute F with the best inliers
    bestF = compute8pointAlgorithm(bestInliers);

    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();

    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(bestInliers[i]);
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const Matrix<float>& F) {
    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;
        // --------------- TODO ------------
        Color color(rand()%256,rand()%256,rand()%256);
        drawCircle(x,y,5, color);

        Vector<float> v(3);
        v[0] = x;
        v[1] = y;
        v[2] = 1;

        int a,b;
        //click on the left
        if(x <= I1.width()){
            v = transpose(F) * v;
            a = (int)(-v[2]/v[1]);
            b = (int)(-(v[0]*I2.width()+v[2])/v[1]);
            drawLine(I1.width(),a,I1.width()+I2.width(),(int)b, color);

        }else{
            // homogeneous cordinates
            v[0] = v[0] - I1.width();
            v = F * v;
            a = (int)(-v[2]/v[1]);
            b = (int)(-(v[0]*I1.width()+v[2])/v[1]);
            drawLine(0,a,I1.width(),b, color);

        }

    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    cout << " matches: " << matches.size() << endl;
    click();

    Matrix<float> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);
    }
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
