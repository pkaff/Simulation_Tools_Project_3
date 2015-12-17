within ;
model woodpecker "woodpecker model for project 3 in simulation tools"

 import Modelica.Utilities.Streams;
  parameter Real mS=3.0*10^(-4) "sleeve mass";
  parameter Real Js=5*10^(-9) "sleeve momentum of inertia";
  parameter Real mB=4.5*10^(-3) "bird mass";
  parameter Real Jb=7*10^(-7) "bird momentum of inertia";
  parameter Real r0=2.5*10^(-3) "bar radius";
  parameter Real rS=3.1*10^(-3) "inner raduis of sleeve";
  parameter Real hS=5.8*10^(-3) "half height of the sleeve";
  parameter Real lS=1*10^(-2)
    "distance between sleeve and spring rotation axis";
  parameter Real lG=1.5*10^(-2) "distance between brid and pring rotation axis";
  parameter Real lB=2.01*10^(-2) "beak x-coordinate in bird's local system";
  parameter Real hB=2.0*10^(-2) "beak y-coordinate in bird's local system";
  parameter Real cP=5.6*10^(-3) "spring konstant";
  parameter Real g=9.81 "gravity";

  Real z(start=0.0) "z-coorcdinate for sleeve";
  Real phiS(start=-0.10344) "sleeve angle";
  Real phiB(start=-0.65) "bird angle";

  Real zdot(start=0.0) "first deriviative of z";
  Real phiSdot(start=0.0) "first deriviative of phiS";
  Real phiBdot(start=0) "first deriviative of phiB";

  Real lambda1(start=-0.6911);
  Real lambda2(start=-0.1416);

  Integer state(start=2, fixed=true) "state of system";
  Integer count4(start=0) "number of times to state 4";

equation
    der(z)=zdot;
    der(phiS)=phiSdot;
    der(phiB)=phiBdot;

    if state==1 then
      //Streams.print("state = " + String(state));
      (mS+mB)*der(zdot)+mB*lS*der(phiSdot)+mB*lG*der(phiBdot) = -(mS+mB)*g;
      (mB*lS)*der(zdot)+(Js+mB*lS^2)*der(phiSdot)+(mB*lS*lG)*der(phiBdot) = cP*(phiB-phiS)-mB*lS*g-lambda1;
      mB*lG*der(zdot)+mB*lS*lG*der(phiSdot)+(Jb+mB*lG^2)*der(phiBdot) = cP*(phiS-phiB)-mB*lG*g-lambda2;
      lambda1=0;
      lambda2=0;
    elseif state==2 then
      //Streams.print("state = " + String(state));
      (mS+mB)*der(zdot)+mB*lS*der(phiSdot)+mB*lG*der(phiBdot) = -(mS+mB)*g-lambda2;
      (mB*lS)*der(zdot)+(Js+mB*lS^2)*der(phiSdot)+(mB*lS*lG)*der(phiBdot) = cP*(phiB-phiS)-mB*lS*g-hS*lambda1-rS*lambda2;
      mB*lG*der(zdot)+mB*lS*lG*der(phiSdot)+(Jb+mB*lG^2)*der(phiBdot) = cP*(phiS-phiB)-mB*lG*g;
      //(rS-r0)+hS*phiS=0;
      0=hS*der(phiSdot);
      der(zdot)+rS*der(phiSdot)=0;
    else
      //Streams.print("state = " + String(state));
      (mS+mB)*der(zdot)+mB*lS*der(phiSdot)+mB*lG*der(phiBdot) = -(mS+mB)*g-lambda2;
      (mB*lS)*der(zdot)+(Js+mB*lS^2)*der(phiSdot)+(mB*lS*lG)*der(phiBdot) = cP*(phiB-phiS)-mB*lS*g+hS*lambda1-rS*lambda2;
      mB*lG*der(zdot)+mB*lS*lG*der(phiSdot)+(Jb+mB*lG^2)*der(phiBdot) = cP*(phiS-phiB)-mB*lG*g;
      //(rS-r0)-hS*phiS=0;
      0=-hS*der(phiSdot);
      der(zdot)+rS*der(phiSdot)=0;
    end if;

algorithm
  //Streams.print("start: " + String(state));
  //Streams.print("lambda1 = " + String(lambda1));
    when state==1 and phiBdot<0 and (hS*phiS+rS-r0)<0 then
      state:=2;
      Streams.print("A");
      reinit(phiBdot, (mB*lG*pre(zdot) + (mB*lS*lG)*pre(phiSdot) + (Jb + mB*lG^2)*pre(phiBdot))/(Jb + mB*lG^2));
      Streams.print("B");
      reinit(zdot, 0);
      Streams.print("C");
      reinit(phiSdot, 0);
      Streams.print("lala: " + String(pre(state))+ " -> " + String(state));
      Streams.print("1 -> 2");
    elsewhen state==1 and phiBdot>0 and (hS*phiS-rS+r0)>0 then
      reinit(phiBdot, (mB*lG*pre(zdot) + (mB*lS*lG)*pre(phiSdot) + (Jb + mB*lG^2)*pre(phiBdot))/(Jb + mB*lG^2));
      reinit(zdot, 0);
      reinit(phiSdot, 0);
      Streams.print(String(state));
      state:=3;
      Streams.print("1 -> 3");
      Streams.print(String(state));
    elsewhen lambda1>10^(-4) then
      if state==2 then
        state:=1;
        Streams.print("2 -> 1 A");
      elseif state==3 and phiBdot<0 then
        state:=1;
        Streams.print("3 -> 1 A");
      end if;
    elsewhen lambda1<-10^(-4) and 1 == 2 then
      if state==2 then
        state:=1;
        Streams.print("2 -> 1 B");
      elseif state==3 and phiBdot<0 then
        state:=1;
        Streams.print("3 -> 1 B");
      end if;
    elsewhen state==3 and phiBdot>0 and hB*phiB-lS-lG+lB+r0>0 then
      reinit(phiBdot,-pre(phiBdot));
      count4 := pre(count4) + 1;
      Streams.print("count4 = " + String(count4));
      Streams.print("3 -> 4 -> 3");
    end when;
    //Streams.print("ldsfsdl");
    //Streams.print("end: " + String(state));

  annotation (uses(Modelica(version="3.2.1")));
end woodpecker;
