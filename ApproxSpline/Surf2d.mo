within ApproxSpline;
package Surf2d
  extends Modelica.Icons.Package;
  function eval
    extends Modelica.Icons.Function;
    input ApproxSpline.Surf2d.Type s "2-D spline object to operate on";
    input Real x "first independent variable where spline is to be evaluated";
    input Real y "second independent variable where spline is to be evaluated";
    output Real z "evaluated value ";
    external "C" z = surf2dEval(s, x, y)
      annotation(Library="mdc_dierckx", __iti_dll="ITI_mdc_dierckx.dll", _iti_dllNoExport=true);
    annotation(derivative = derivative);
  end eval;

  function derivative
    extends Modelica.Icons.Function;
    input ApproxSpline.Surf2d.Type s "2-D spline object to operate on";
    input Real x "first independent variable where spline is to be evaluated";
    input Real y "second independent variable where spline is to be evaluated";
    input Real ddx "value of the 1st partial derivative";
    input Real ddy "value of the 2nd partial derivative";
    output Real dz "evaluated value ";
    external "C" dz = surf2dDer(s, x, y, ddx, ddy)
      annotation(Library="mdc_dierckx", __iti_dll="ITI_mdc_dierckx.dll", _iti_dllNoExport=true);
  end derivative;

  function getNumberOfKnotsX
    extends Modelica.Icons.Function;
    input ApproxSpline.Surf2d.Type s "2-D spline object to operate on";
    output Integer nx "number of knots in X-dir";
    external "C" nx = surf2dGetNumberOfKnotsX(s)
      annotation(Library="mdc_dierckx", __iti_dll="ITI_mdc_dierckx.dll", _iti_dllNoExport=true);
  end getNumberOfKnotsX;

  function getNumberOfKnotsY
    extends Modelica.Icons.Function;
    input ApproxSpline.Surf2d.Type s "2-D spline object to operate on";
    output Integer ny "number of knots in Y-dir";
    external "C" ny = surf2dGetNumberOfKnotsY(s)
      annotation(Library="mdc_dierckx", __iti_dll="ITI_mdc_dierckx.dll", _iti_dllNoExport=true);
  end getNumberOfKnotsY;

class Type "External object of a 2-D (X,Y) -> Z approximating spline surface"
  extends ExternalObject;
  extends Icons.External;

  function constructor "create a 2D-Spline object"
    extends Modelica.Icons.Function;
    input Boolean rectangular =  false
        "set true if input data is on rectangular grid, otherwise assume scattered data";
    input Real data[:,:] "data to be fitted (see doku)";
    input Real s(min=0) = 0 "smoothing factor (s>0)";
    input Integer kx(min=1,max=5) = 3 "degree of X-dir spline (should be odd)";
    input Integer ky(min=1,max=5) = 3 "degree of Y-dir spline (should be odd)";
    input Real x_lim[2] = {
      if rectangular then  ApproxSpline.Utilities.vmin(data[2:size(data,1),1]) else  ApproxSpline.Utilities.vmin(data[:,1]),
      if rectangular then  ApproxSpline.Utilities.vmax(data[2:size(data,1),1]) else  ApproxSpline.Utilities.vmax(data[:,1])}
        "X-dir boundaries of interpolation interval";
    input Real y_lim[2] = {
      if rectangular then  ApproxSpline.Utilities.vmin(data[1,2:size(data,2)]) else  ApproxSpline.Utilities.vmin(data[:,2]),
      if rectangular then  ApproxSpline.Utilities.vmax(data[1,2:size(data,2)]) else  ApproxSpline.Utilities.vmax(data[:,2])}
        "Y-dir boundaries of interpolation interval";
    input Real tx[:] = fill(0,0)
        "array of X-dir knot positions (if given together with ty, a least square spline is generated)";
    input Real ty[:] = fill(0,0)
        "array of Y-dir knot positions (if given together with ty, a least square spline is generated)";
    output Type surf "2-D spline surface data structure";
    external "C" surf = surf2dNew(rectangular, data, size(data,1), size(data,2), x_lim, y_lim, kx, ky, s, tx, ty, size(tx,1), size(ty,1))
      annotation(Library="mdc_dierckx", __iti_dll="ITI_mdc_dierckx.dll", _iti_dllNoExport=true);
    annotation (Documentation(info="<html>
<h4>Generate smooth 2D spline surface</h4>
This function generates a 2D-surface object by approximating spline fitting.
Input data may be scattered or on a rectangular grid.

<h5>Rectangular data</h5>
In case of rectangular data points, set boolean parameter isRect to true and provide the data as array in the common Modelica table form:

<pre>

  0 | y1 | y2 | y3 | ...
  --+----+----+----+-----
 x1 | z11| z12| z13| ...
 x2 | z21| z22| z23| ...
 x3 | z31| z32| z33| ...
 ...| ...| ...| ...| ...

</pre>

<h5>Scattered data</h5>
In case of scattered data points, set parameter isRect to false and provide the data as 3- or 4-column array:

<pre>

 x1 | y1 | z1 | w1 |
 x2 | y2 | z2 | w2 |
 x3 | y3 | z3 | w3 |
 ...| ...| ...| .. |

or

 x1 | y1 | z1 |
 x2 | y2 | z2 |
 x3 | y3 | z3 |
 ...| ...| ...|

</pre>



The last column w[i] is the optional weight of the data points. If no weights are given, any data points are weighted equally.

</html>"));
  end constructor;

  function destructor "remove a 2D-Spline object"
    extends Modelica.Icons.Function;
    input Type surf "2-D spline surface data structure";
    external "C" surf2dDel(surf)
      annotation(Library="mdc_dierckx", __iti_dll="ITI_mdc_dierckx.dll", _iti_dllNoExport=true);
  end destructor;

    annotation (Documentation(info="

"));
end Type;
end Surf2d;
