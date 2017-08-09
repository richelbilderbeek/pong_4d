/*


Pin |Connect to 
----|-----
A0  |Joystick 1 X
A1  |Joystick 1 Y
A2  |Joystick 1 switch
A3  |Joystick 2 X
A4  |Joystick 2 Y
A5  |Joystick 2 switch

*/

void setup() 
{
  pinMode(A0, INPUT);
  pinMode(A1, INPUT);
  pinMode(A2, INPUT);
  pinMode(A3, INPUT);
  pinMode(A4, INPUT);
  pinMode(A5, INPUT);
  Serial.begin(9600);
}

void loop() 
{
  Serial.write(0);
  Serial.write(map(analogRead(A0), 0, 1023, 1, 255));
  Serial.write(map(analogRead(A1), 0, 1023, 1, 255));
  Serial.write(map(analogRead(A2), 0, 1023, 1, 255));
  Serial.write(map(analogRead(A3), 0, 1023, 1, 255));
  Serial.write(map(analogRead(A4), 0, 1023, 1, 255));
  Serial.write(map(analogRead(A5), 0, 1023, 1, 255));
  delay(20);
}

