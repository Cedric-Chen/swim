advance fluid
=============

tn -> tn+1
+ input 
    + u\_fluid, u\_beam, char fn @tn
+ output 
    + vorticity @tn+1


+ u\_fluid @tn, u\_beam @tn -> u\_lamda
+ u\_lamda -> vorticity
+ vorticity @tn -> vorticity @tn+1

advance solid
=============

tn -> tn+1
+ input
    + vorticity @tn+1
    + angle, velangle @tn-1, tn
+ output
    + angle @tn+1, velangle @tn+1
    + char fn @tn+1, u\_fluid, u\_beam @tn+1


+ kill vorticity @tn+1 at the boundary 
+ vorticity -> u\_fluid @tn+1
+ u\_fluid @tn+1 -> pressure @tn+1

+ use secant method:
    angle, velangle @tn-1, @tn -> angle, velangle @tn+1
    + estimate value for angle, velangle @tn+1;
    + beam.update(), beam.update\_pressure -> pressureOnBody @tn+1
    + solve for T @tn+1
    + solve for g, g is a funtion of angle, velangle @tn-1, tn, tn+1


class advanceBeam
-----------------
(currPressure, oldAngle, oldVelAngle, currAngle, currVelAngle, currT,
    newAngle, newVelAngle, newT)

+ member variable
    + currPressure
    + oldAngle, oldVelAngle, currAngle, currVelAngle, currT
    + newAngle, newVelAngle, newT
+ member function
    + operator()
        + the main flow of the program
        + broyden(computeF, currangel, newangle, newvelangle)
    + computeF(estimatedAngle, estimatedVelAngle)
        + beam(currangle, currvelangle)
        + pressureOnBody = vector{..}
        + beam.update_pressure(pressureOnBody)
        + broygen(computeG, estimatedT, tempT)
        + return ----
    + computeG(estimatedT)
        + return ----
    + broyden(estimated, result)
+ __mistake__
    + I think I should combine class advanceBeam and class beam together.
        + they two are both properties of beam.
        + advance need to use member function of class beam: updatePressure.

Class Beam
==========

member variable
---------------
+ input
    + goemetry: width, length, locHead, velHead
    + computational setting: epsilon, numofElem
    + fluid knowledge: pressure
    + dynamic knowledge: oldAngle, oldVelAngle, currAngle, currVelAngle, currT
+ internal variable
    + deltaFunc
    + pressureOnBody
+ output
    + charFunc
    + velBeam
    + newangle, newvelangle, newT

member function
---------------
+ public
    + constructor(geometry and computational setting)
    + advance(pressure, 3 angle, 3 velAngle, 2 T, velBeam, charFunc)
        + Broyden()
+ private
    + update(angle, velangle)
        + ---
        + update_pressure
