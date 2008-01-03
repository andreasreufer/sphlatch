import processing.opengl.*;
import processing.pdf.*;

float boxSize = 40;
float margin = boxSize*2;
float depth = 500;
color boxFill;

String[] inputStrings;
int noNodes = 0;

boolean[] isCell;
float[] xPos, yPos, zPos, cellSize;

void setup(){
  inputStrings = loadStrings("treedump20.txt");
  //inputStrings = loadStrings("treedump50.txt");
  
  noNodes = inputStrings.length;
  isCell = new boolean[noNodes];
  xPos = new float[noNodes];
  yPos = new float[noNodes];
  zPos = new float[noNodes];
  cellSize = new float[noNodes];
  
  for (int i = 0; i < noNodes; i++) {
    String[] comp = split(inputStrings[i], "  ");
    String[] attrs = split(comp[0], "");
    
    xPos[i] = float(comp[2]);
    yPos[i] = float(comp[3]);
    zPos[i] = float(comp[4]);
  
    String[] pmatch = match( attrs[1], "P" );  
    if ( pmatch == null ) {
      isCell[i] = true;
      cellSize[i] = float(comp[5]);
    } else {
      isCell[i] = false;
      cellSize[i] = 0.;
    }    
  }
  
  size(600, 600, PDF, "output.pdf");
  //size(600, 600, P2D);
  noStroke();
}

void draw(){
 
  // center and spin grid
  
  /*translate(width/2, height/2);
  rotateY(frameCount*PI/300);
  rotateX(frameCount*PI/300);
  translate(-depth/2, -depth/2, -depth/2);*/
  
  //background(0);
  background(255);
 
  /*strokeWeight(2);
  translate(-0.1*depth, -0.1*depth, -0.1*depth);
  stroke(#FF0000, 200);
  line(0,0,0,1*depth,0,0);
  stroke(#00FF00, 200);
  line(0,0,0,0,1*depth,0);
  stroke(#0000FF, 200);
  line(0,0,0,0,0,1*depth);
  translate(0.1*depth, 0.1*depth, 0.1*depth);*/
  translate(50, 50);    
  
  for (int i = 0; i < noNodes; i++) {
    translate(depth*xPos[i], depth*yPos[i]);
    
    if ( isCell[i] ) {

      stroke(#000000, 255);
      noFill();
      //fill(#00FF00, 5);
      //box( depth*cellSize[i] );
      float dx = 0.5*depth*cellSize[i];
      quad( -dx, dx,
          dx, dx,
          dx, -dx,
          -dx, -dx );
    } else {
      stroke(#00FF00, 255);
      fill(#00FF00, 255);
      ellipse(0,0,5,5);
      //point(0,0);
      //sphereDetail(1);
      //sphere(1);
    }
     
    translate(-depth*xPos[i], -depth*yPos[i]);   
  }
  
  exit();
}

