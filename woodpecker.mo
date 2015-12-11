within ;
model woodpecker "woodpecker model for project 3 in simulation tools"

 import Modelica.Utilities.Streams;
  parameter Real mS=3.0*10^4 "sleeve mass";
  parameter Real Js=5*10^(-9) "sleeve momentum of inertia";
  parameter Real mB=4.5*10^(-3) "bird mass";
  parameter Real Jb=7*10^(-7) "bird momentum of inertia";
  parameter Real r0=2.5*10^(-3) "bar radius";
  parameter Real rS=3.1*10^(-3) "inner raduis of sleeve";
  parameter Real hS=2* 10^(-2) "half height of the sleeve";
  parameter Real lS=1*10^(-2)
    "distance between sleeve and spring rotation axis";
  parameter Real lG=1.5*10^(-2) "distance between brid and pring rotation axis";
  parameter Real lB=2.01*10^(-2) "beak x-coordinate in bird's local system";
  parameter Real hB=2.0*10^(-2) "beak y-coordinate in bird's local system";
  parameter Real cP=5.6*10^(-3) "spring konstant";
  parameter Real g=9.81 "gravity";

  Real z(start=10) "z-coorcdinate for sleeve";
  Real phiS(start=0.08) "sleeve angle";
  Real phiB(start=-1) "bird angle";

  Real zdot "first deriviative of z";
  Real phiSdot "first deriviative of phiS";
  Real phiBdot "first deriviative of phiB";

  Real lambda1;
  Real lambda2;

  Integer state(start=1) "state of system";

equation
    der(z)=zdot;
    der(phiS)=phiSdot;
    der(phiB)=phiBdot;

    if state==1 then
      (mS+mB)*der(zdot)+mB*lS*der(phiSdot)+mB*lG*der(phiBdot) = -(mS+mB)*g;
      (mB*lS)*der(zdot)+(Js+mB*lS^2)*der(phiSdot)+(mB*lS*lG)*der(phiBdot) = cP*(phiB-phiS)-mB*lS*g-lambda1;
      mB*lG*der(zdot)+mB*lS*lG*der(phiSdot)+(Jb+mB*lG^2)*der(phiBdot) = cP*(phiS-phiB)-mB*lG*g-lambda2;
      der(lambda1)=0;
      der(lambda2)=0;
    elseif state==2 then
      (mS+mB)*der(zdot)+mB*lS*der(phiSdot)+mB*lG*der(phiBdot) = -(mS+mB)*g-lambda2;
      (mB*lS)*der(zdot)+(Js+mB*lS^2)*der(phiSdot)+(mB*lS*lG)*der(phiBdot) = cP*(phiB-phiS)-mB*lS*g-hS*lambda1-rS*lambda2;
      mB*lG*der(zdot)+mB*lS*lG*der(phiSdot)+(Jb+mB*lG^2)*der(phiBdot) = cP*(phiS-phiB)-mB*lG*g;
      //(rS-r0)+hS*phiS=0;
      0=hS*der(phiSdot)*phiSdot;
      zdot+rS*phiSdot=0;
    else
      (mS+mB)*der(zdot)+mB*lS*der(phiSdot)+mB*lG*der(phiBdot) = -(mS+mB)*g-lambda2;
      (mB*lS)*der(zdot)+(Js+mB*lS^2)*der(phiSdot)+(mB*lS*lG)*der(phiBdot) = cP*(phiB-phiS)-mB*lS*g+hS*lambda1-rS*lambda2;
      mB*lG*der(zdot)+mB*lS*lG*der(phiSdot)+(Jb+mB*lG^2)*der(phiBdot) = cP*(phiS-phiB)-mB*lG*g;
      //(rS-r0)-hS*phiS=0;
      0=-hS*der(phiSdot)*phiSdot;
      zdot+rS*phiSdot=0;
    end if;

algorithm

    if state==1 and phiBdot<0 and (hS*phiS+rS-r0)<0 then
      state:=2;
      Streams.print("1 -> 2");
    elseif state==1 and phiBdot>0 and (hS*phiS-rS+r0)>0 then
      state:=3;
      Streams.print("1 -> 3");
    end if;

       when lambda1>0 then
          if state==2 then
          state:=1;
          Streams.print("2 -> 1");
          end if;
       end when;
      when lambda1<0 then
        if state==2 then
          state:=1;
          Streams.print("2 -> 1");
        end if;
      end when;

       when lambda1>0 then
        if state==3 and phiBdot<0 then
          state:=1;
          Streams.print("3 -> 1");
        end if;
      end when;

      when lambda1<0 then
        if state==3 and phiBdot<0 then
        state:=1;
        Streams.print("3 -> 1");
        end if;
      end when;

    if state==3 and phiBdot>0 and hB*phiB-lS-lG+lB+r0>0 then
      state:=4;
      Streams.print("3 -> 4");
    end if;

    when state==4 then
      state:=3;
      reinit(phiBdot,-pre(phiBdot));
      Streams.print("4 -> 3");
    end when;

  annotation (uses(Modelica(version="3.2.1")));
end woodpecker;
