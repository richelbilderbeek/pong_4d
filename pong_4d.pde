import processing.serial.*;

import ddf.minim.*;
AudioPlayer player;
Minim minim;//audio context

//import processing.sound.*;
//SoundFile file;

// Arena is 500 x 500 x 500
// Center of the arena is (0, 0, 0)
final float speed_z = 0.01;
//final float speed_z = 0.001;

final float player1_r = 50;
final float player2_r = 50;
final float bal_r = 50.0;
final float minx = -500.0;
final float maxx =  500.0;
final float miny = -250.0;
final float maxy =  250.0;

final float bal_x = 0.0;
final float bal_y = 0.0;
final float bal_dx = 5.0;
final float bal_dy = 0.1;
final float player1_x = minx + (1.0 * player1_r);
final float player1_y = 0.0;
final float player1_dy = 0.0;
final float player2_x = maxx - (1.0 * player2_r);
final float player2_y = 0.0;
final float player2_dy = 0.0;
final float text_size = 128;
float hoek_z = 0.0;
int score1 = 0;
int score2 = 0;

Ball player_1 = new Ball(
  new PVector(player1_x, player1_y), 
  new PVector(0.0, 0.0), 
  player1_r);

Ball player_2 = new Ball(
  new PVector(player2_x, player2_y), 
  new PVector(0.0, 0.0), 
  player2_r);

Ball bal = new Ball(
  new PVector(bal_x, bal_y), 
  new PVector(bal_dx, bal_dy), 
  bal_r);

Serial serial_port;



/// standalone means without Arduino
final boolean standalone = true;

void setup()
{
  fullScreen(P3D);
  textSize(text_size);
  if (!standalone)
  {
    serial_port = new Serial(this, Serial.list()[0], 9600);
  }

  minim = new Minim(this);
  player = minim.loadFile("EddiesTwister.mp3");
  player.play();
}

void draw()
{
  if (!standalone)
  {
    while (serial_port.available() > 0) 
    {
      final int spacer_value = serial_port.read();
      if (spacer_value != 0) continue; 
      if (serial_port.available() == 0) break;
      player_1.set_y(map(serial_port.read(), 1, 255, miny, maxy));
      if (serial_port.available() == 0) break;
      player_2.set_y(map(serial_port.read(), 1, 255, miny, maxy));
    }
  }
  background(0);
  noStroke();
  directionalLight(255, 0, 0, -cos(hoek_z * 5.0), -sin(hoek_z * 5.0), 0);
  directionalLight(0, 255, 0, 0, -cos(hoek_z * 7.0), -sin(hoek_z * 7.0));
  directionalLight(0, 0, 255, -sin(hoek_z * 3.0), 0, -cos(hoek_z * 3.0));

  translate(width / 2, height / 2);
  rotateX(hoek_z * 0.2);
  rotateY(hoek_z * 0.3);
  rotateZ(hoek_z);

  player_1.move();
  player_2.move();
  bal.move();

  bal.check_collision(player_1);
  bal.check_collision(player_2);

  player_1.display();
  player_2.display();
  bal.display();

  draw_surrounding();

  fill(255, 128, 128);
  text(score1, -text_size, text_size / 2);

  fill(128, 128, 255);
  text(score2,  text_size, text_size / 2);

  hoek_z += speed_z;
}

void draw_surrounding()
{
  //Draw surrounding
  for (float i = 0; i < 100.0; i += 1.0)
  {
    final float block_x = minx + ((i / 100.0) * (maxx - minx));
    final float block_y = miny ;
    translate(block_x, block_y);
    box(10);
    translate(-block_x, -block_y);
  }
  for (float i = 0; i < 100.0; i += 1.0)
  {
    final float block_x = minx + ((i / 100.0) * (maxx - minx));
    final float block_y = maxy;
    translate(block_x, block_y);
    box(10);
    translate(-block_x, -block_y);
  }
  for (float i = 0; i < 100.0; i += 1.0)
  {
    final float block_x = minx;
    final float block_y = miny + ((i / 100.0) * (maxy - miny));
    translate(block_x, block_y);
    box(10);
    translate(-block_x, -block_y);
  }
  for (float i = 0; i < 100.0; i += 1.0)
  {
    final float block_x = maxx;
    final float block_y = miny + ((i / 100.0) * (maxy - miny));
    translate(block_x, block_y);
    box(10);
    translate(-block_x, -block_y);
  }
}

float get_distance(float dx, float dy)
{
  return sqrt((dx * dx) + (dy * dy));
}

void keyPressed() {

  if (key == CODED) 
  {
    if (keyCode == UP) 
    {
      player_2.change_dy(-1);
    } 
    if (keyCode == DOWN) 
    {
      player_2.change_dy(1);
    }
  }
  if (key == 'w') 
  {
    player_1.change_dy(-1);
  }
  if (key == 's') 
  {
    player_1.change_dy( 1);
  }
}


// Use Ball from https://processing.org/examples/circlecollision.html
class Ball 
{
  PVector m_position;
  PVector m_velocity;
  float m_radius;

  Ball(PVector position, PVector velocity, float radius) 
  {
    m_position = position;
    m_velocity = velocity;
    m_radius = radius;
  }

  final float get_radius() { 
    return m_radius;
  }

  void set_velocity_to_zero()
  {
    m_velocity.x = 0.0; 
    m_velocity.y = 0.0;
  }
  void set_y(float y) 
  {
    m_position.y = y;
  }
  void change_dy(float change)
  {
    m_velocity.y += change;
  }
  void move() 
  {
    m_position.add(m_velocity);

    if (m_position.x > maxx - m_radius) 
    {
      m_velocity.x = -abs(m_velocity.x);
      ++score1;
    } 
    else if (m_position.x < minx + m_radius) 
    {
      m_velocity.x = abs(m_velocity.x);
      ++score2;
    } 
    else if (m_position.y > maxy - m_radius) 
    {
      m_velocity.y = -abs(m_velocity.y);
    } 
    else if (m_position.y < miny + m_radius) 
    {
      m_velocity.y = abs(m_velocity.y);
    }
  }

  void check_collision(final Ball other) 
  {

    // Get distances between the balls components
    final PVector distanceVect = PVector.sub(other.m_position, m_position);

    // Calculate magnitude of the vector separating the balls
    final float distanceVectMag = distanceVect.mag();

    // Minimum distance before they are touching
    final float minDistance = m_radius + other.get_radius();

    //No collission
    if (distanceVectMag > minDistance) return;

    final float distanceCorrection = (minDistance-distanceVectMag)/2.0;
    final PVector d = distanceVect.copy();
    final PVector correctionVector = d.normalize().mult(distanceCorrection);
    
    //RJCB: I don't think I need this
    //other.m_position.add(correctionVector);
    m_position.sub(correctionVector);

    // get angle of distanceVect
    final float theta  = distanceVect.heading();
    // precalculate trig values
    final float sine = sin(theta);
    final float cosine = cos(theta);

    /* bTemp will hold rotated ball positions. You 
     just need to worry about bTemp[1] position*/
    PVector[] bTemp = {
      new PVector(), new PVector()
    };

    /* this ball's position is relative to the other
     so you can use the vector between them (bVect) as the 
     reference point in the rotation expressions.
     bTemp[0].position.x and bTemp[0].position.y will initialize
     automatically to 0.0, which is what you want
     since b[1] will rotate around b[0] */
    bTemp[1].x  = cosine * distanceVect.x + sine * distanceVect.y;
    bTemp[1].y  = cosine * distanceVect.y - sine * distanceVect.x;

    // rotate Temporary velocities
    PVector[] vTemp = {
      new PVector(), new PVector()
    };

    vTemp[0].x  = cosine * m_velocity.x + sine * m_velocity.y;
    vTemp[0].y  = cosine * m_velocity.y - sine * m_velocity.x;
    
    //Instead of using the player's velocity, use the opposite velocity of the ball 
    //vTemp[1].x  = cosine * other.m_velocity.x + sine * other.m_velocity.y;
    //vTemp[1].y  = cosine * other.m_velocity.y - sine * other.m_velocity.x;
    vTemp[1].x  = cosine * -m_velocity.x + sine * -m_velocity.y;
    vTemp[1].y  = cosine * -m_velocity.y - sine * -m_velocity.x;

    /* Now that velocities are rotated, you can use 1D
     conservation of momentum equations to calculate 
     the final velocity along the x-axis. */
    PVector[] vFinal = {  
      new PVector(), new PVector()
    };

    // RJCB: No idea what these m's are for
    //final float m = m_radius;
    //final float other_m = other.get_radius();
    //final float m = m_radius * 0.1;
    //final float other_m = other.get_radius() * 0.1;

    // final rotated velocity for b[0]

    vFinal[0].x = ((m_radius - other.get_radius()) * vTemp[0].x + 2 * other.get_radius() * vTemp[1].x) / (m_radius + other.get_radius());
    vFinal[0].y = vTemp[0].y;

    // final rotated velocity for b[0]
    vFinal[1].x = ((other.get_radius() - m_radius) * vTemp[1].x + 2 * m_radius * vTemp[0].x) / (m_radius + other.get_radius());
    vFinal[1].y = vTemp[1].y;

    // hack to avoid clumping
    bTemp[0].x += vFinal[0].x;
    bTemp[1].x += vFinal[1].x;

    /* Rotate ball positions and velocities back
     Reverse signs in trig expressions to rotate 
     in the opposite direction */
    // rotate balls
    PVector[] bFinal = { 
      new PVector(), new PVector()
    };

    bFinal[0].x = cosine * bTemp[0].x - sine * bTemp[0].y;
    bFinal[0].y = cosine * bTemp[0].y + sine * bTemp[0].x;
    bFinal[1].x = cosine * bTemp[1].x - sine * bTemp[1].y;
    bFinal[1].y = cosine * bTemp[1].y + sine * bTemp[1].x;

    //Nope, other players do not move
    // update balls to screen position
    //other.m_position.x = m_position.x + bFinal[1].x;
    //other.m_position.y = m_position.y + bFinal[1].y;

    m_position.add(bFinal[0]);

    // update velocities
    m_velocity.x = cosine * vFinal[0].x - sine * vFinal[0].y;
    m_velocity.y = cosine * vFinal[0].y + sine * vFinal[0].x;
    
    //Nope, other players do not move
    //other.m_velocity.x = cosine * vFinal[1].x - sine * vFinal[1].y;
    //other.m_velocity.y = cosine * vFinal[1].y + sine * vFinal[1].x;

  }

  void display() 
  {
    translate(m_position.x, m_position.y);
    sphere(m_radius);
    translate(-m_position.x, -m_position.y);
  }
}




// From http://www.plasmaphysics.org.uk/programs/coll3d_cpp.htm:
/*****************************************************************************
//   This program is a 'remote' 3D-collision detector for two balls on linear
//   trajectories and returns, if applicable, the location of the collision for 
//   both balls as well as the new velocity vectors (assuming a partially elastic
//   collision as defined by the restitution coefficient).
//
//   All variables apart from 'error' are of Double Precision Floating Point type.
//
//   The Parameters are:
//
//    R    (restitution coefficient)  between 0 and 1 (1=perfectly elastic collision)
//    m1    (mass of ball 1)
//    m2    (mass of ball 2)
//    r1    (radius of ball 1)
//    r2    (radius of ball 2)
//  & x1    (x-coordinate of ball 1) 
//  & y1    (y-coordinate of ball 1)          
//  & z1    (z-coordinate of ball 1) 
//  & x2    (x-coordinate of ball 2)              
//  & y2    (y-coordinate of ball 2)         
//  & z2    (z-coordinate of ball 2)         
//  & vx1   (velocity x-component of ball 1) 
//  & vy1   (velocity y-component of ball 1)
//  & vz1   (velocity z-component of ball 1)          
//  & vx2   (velocity x-component of ball 2)         
//  & vy2   (velocity y-component of ball 2)
//  & vz2   (velocity z-component of ball 2)
//  & error (int)     (0: no error 
//                     1: balls do not collide
//                     2: initial positions impossible (balls overlap))
//
//   Note that the parameters with an ampersand (&) are passed by reference,
//   i.e. the corresponding arguments in the calling program will be updated 
//   (the positions and velocities however only if 'error'=0).
//   All variables should have the same data types in the calling program
//   and all should be initialized before calling the function.
//
//   This program is free to use for everybody. However, you use it at your own
//   risk and I do not accept any liability resulting from incorrect behaviour.
//   I have tested the program for numerous cases and I could not see anything 
//   wrong with it but I can not guarantee that it is bug-free under any 
//   circumstances.
//
//   I would appreciate if you could report any problems to me
//   (for contact details see  http://www.plasmaphysics.org.uk/feedback.htm ).
//
//   Thomas Smid   February 2004
//                 December 2005 (a few minor changes to improve speed)
//                 December 2009 (generalization to partially inelastic collisions)
//                 July     2011 (fix for possible rounding errors)
//******************************************************************************

   
    void collision3D(double R, double m1, double m2, double r1, double r2,
                     double& x1, double& y1,double& z1,
                     double& x2, double& y2, double& z2,
                     double& vx1, double& vy1, double& vz1,
                     double& vx2, double& vy2, double& vz2,
                     int& error)     {


       double  pi,r12,m21,d,v,theta2,phi2,st,ct,sp,cp,vx1r,vy1r,vz1r,fvz1r,
             thetav,phiv,dr,alpha,beta,sbeta,cbeta,dc,sqs,t,a,dvz2,
         vx2r,vy2r,vz2r,x21,y21,z21,vx21,vy21,vz21,vx_cm,vy_cm,vz_cm;

//     **** initialize some variables ****
       pi=acos(-1.0E0);
       error=0;
       r12=r1+r2;
       m21=m2/m1;
       x21=x2-x1;
       y21=y2-y1;
       z21=z2-z1;
       vx21=vx2-vx1;
       vy21=vy2-vy1;
       vz21=vz2-vz1;
       
       vx_cm = (m1*vx1+m2*vx2)/(m1+m2) ;
       vy_cm = (m1*vy1+m2*vy2)/(m1+m2) ;
       vz_cm = (m1*vz1+m2*vz2)/(m1+m2) ;  

     
//     **** calculate relative distance and relative speed ***
       d=sqrt(x21*x21 +y21*y21 +z21*z21);
       v=sqrt(vx21*vx21 +vy21*vy21 +vz21*vz21);
       
//     **** return if distance between balls smaller than sum of radii ****
       if (d<r12) {error=2; return;}
       
//     **** return if relative speed = 0 ****
       if (v==0) {error=1; return;}
       

//     **** shift coordinate system so that ball 1 is at the origin ***
       x2=x21;
       y2=y21;
       z2=z21;
       
//     **** boost coordinate system so that ball 2 is resting ***
       vx1=-vx21;
       vy1=-vy21;
       vz1=-vz21;

//     **** find the polar coordinates of the location of ball 2 ***
       theta2=acos(z2/d);
       if (x2==0 && y2==0) {phi2=0;} else {phi2=atan2(y2,x2);}
       st=sin(theta2);
       ct=cos(theta2);
       sp=sin(phi2);
       cp=cos(phi2);


//     **** express the velocity vector of ball 1 in a rotated coordinate
//          system where ball 2 lies on the z-axis ******
       vx1r=ct*cp*vx1+ct*sp*vy1-st*vz1;
       vy1r=cp*vy1-sp*vx1;
       vz1r=st*cp*vx1+st*sp*vy1+ct*vz1;
       fvz1r = vz1r/v ;
       if (fvz1r>1) {fvz1r=1;}   // fix for possible rounding errors
          else if (fvz1r<-1) {fvz1r=-1;} 
       thetav=acos(fvz1r);
       if (vx1r==0 && vy1r==0) {phiv=0;} else {phiv=atan2(vy1r,vx1r);}

                    
//     **** calculate the normalized impact parameter ***
       dr=d*sin(thetav)/r12;


//     **** return old positions and velocities if balls do not collide ***
       if (thetav>pi/2 || fabs(dr)>1) {
           x2=x2+x1;
           y2=y2+y1;
           z2=z2+z1;
           vx1=vx1+vx2;
           vy1=vy1+vy2;
           vz1=vz1+vz2;
           error=1;
           return;
        }
       
//     **** calculate impact angles if balls do collide ***
       alpha=asin(-dr);
       beta=phiv;
       sbeta=sin(beta);
       cbeta=cos(beta);
        
       
//     **** calculate time to collision ***
       t=(d*cos(thetav) -r12*sqrt(1-dr*dr))/v;

     
//     **** update positions and reverse the coordinate shift ***
       x2=x2+vx2*t +x1;
       y2=y2+vy2*t +y1;
       z2=z2+vz2*t +z1;
       x1=(vx1+vx2)*t +x1;
       y1=(vy1+vy2)*t +y1;
       z1=(vz1+vz2)*t +z1;
        
 
       
//  ***  update velocities ***

       a=tan(thetav+alpha);

       dvz2=2*(vz1r+a*(cbeta*vx1r+sbeta*vy1r))/((1+a*a)*(1+m21));
       
       vz2r=dvz2;
       vx2r=a*cbeta*dvz2;
       vy2r=a*sbeta*dvz2;
       vz1r=vz1r-m21*vz2r;
       vx1r=vx1r-m21*vx2r;
       vy1r=vy1r-m21*vy2r;

       
//     **** rotate the velocity vectors back and add the initial velocity
//           vector of ball 2 to retrieve the original coordinate system ****
                     
       vx1=ct*cp*vx1r-sp*vy1r+st*cp*vz1r +vx2;
       vy1=ct*sp*vx1r+cp*vy1r+st*sp*vz1r +vy2;
       vz1=ct*vz1r-st*vx1r               +vz2;
       vx2=ct*cp*vx2r-sp*vy2r+st*cp*vz2r +vx2;
       vy2=ct*sp*vx2r+cp*vy2r+st*sp*vz2r +vy2;
       vz2=ct*vz2r-st*vx2r               +vz2;
        

//     ***  velocity correction for inelastic collisions ***

       vx1=(vx1-vx_cm)*R + vx_cm;
       vy1=(vy1-vy_cm)*R + vy_cm;
       vz1=(vz1-vz_cm)*R + vz_cm;
       vx2=(vx2-vx_cm)*R + vx_cm;
       vy2=(vy2-vy_cm)*R + vy_cm;
       vz2=(vz2-vz_cm)*R + vz_cm;  

       return;
}

*/