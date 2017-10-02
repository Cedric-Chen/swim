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
