within ApproxSpline;
package Blocks "Sublibrary of ready to use blocks"

  block Curve1d "1-D spline approximation X->Y"
    extends Modelica.Blocks.Interfaces.SISO;
    extends Icons.Curve1d;

    parameter Real data[:,:]
      "data to be fitted: [x,y,w] (w is optional weight)";
    parameter Real s(min=0) = 0 "smoothing factor (s>=0)";
    parameter Integer k(min=1,max=5) = 3 "degree of spline (should be odd)";
    parameter Boolean periodic = false
      "set true to generate periodic spline curve";
    parameter Real x_lim[2] = { ApproxSpline.Utilities.vmin(data[:,1]), ApproxSpline.Utilities.vmax(data[:,1])}
      "boundaries of interpolation intervall";
    parameter Real t[:] = fill(0,0)
      "array of knot positions (if given, a least square spline is generated)";

  protected
    ApproxSpline.Curve1d.Type spl = ApproxSpline.Curve1d.Type(data=data, s=s, k=k, x_lim=x_lim, t=t);

  equation
    y = ApproxSpline.Curve1d.eval(spl, u);

    annotation (Documentation(info="<html>


<h2>Provide smooth 1D spline curve and interpolate input data</h2>
This blocks generates a 1D-curve object by approximating spline fitting. 
 
 
<h3>Input data</h3>
The data has to be provided as 2- or 3-column array:
 
<pre>
 
 x1 | y1 | w1 |
 x2 | y2 | w2 |
 x3 | y3 | w3 |
 ...| ...| .. |
 
or
 
 x1 | y1 | 
 x2 | y2 | 
 x3 | y3 | 
 ...| ...| 
 
</pre> 
 
 
 
The last column w[i] is the optional weight of the data points. If no weights are given, any data points are weighted equally. 
 

</html>"));
  end Curve1d;

  block Surf2d "2-D spline approximation (X,Y)->Z"
    extends Modelica.Blocks.Interfaces.SI2SO;
    extends Icons.Surf2d;

    parameter Boolean rectangular = false
      "set true if input data is on rectangular grid, otherwise assume scattered data";
    parameter Real data[:,:] "data to be fitted: (see documentation)";
    parameter Real s(min=0) = 0 "smoothing factor (s>=0)";
    parameter Integer kx(min=1,max=5) = 3
      "degree of X-dir spline (should be odd)";
    parameter Integer ky(min=1,max=5) = 3
      "degree of Y-dir spline (should be odd)";
    parameter Real x_lim[2] = {
      if rectangular then  ApproxSpline.Utilities.vmin(data[2:size(data,1),1]) else  ApproxSpline.Utilities.vmin(data[:,1]),
      if rectangular then  ApproxSpline.Utilities.vmax(data[2:size(data,1),1]) else  ApproxSpline.Utilities.vmax(data[:,1])}
      "X-dir boundaries of interpolation intervall";
    parameter Real y_lim[2] = {
      if rectangular then  ApproxSpline.Utilities.vmin(data[1,2:size(data,2)]) else  ApproxSpline.Utilities.vmin(data[:,2]),
      if rectangular then  ApproxSpline.Utilities.vmax(data[1,2:size(data,2)]) else  ApproxSpline.Utilities.vmax(data[:,2])}
      "Y-dir boundaries of interpolation intervall";
    parameter Real tx[:] = fill(0,0)
      "array of X-dir knot positions (if given together with ty, a least square spline is generated)";
    parameter Real ty[:] = fill(0,0)
      "array of Y-dir knot positions (if given together with ty, a least square spline is generated)";

  protected
    ApproxSpline.Surf2d.Type spl = ApproxSpline.Surf2d.Type(rectangular, data, s, kx, ky, x_lim, y_lim, tx, ty);

  equation
    y = ApproxSpline.Surf2d.eval(spl, u1, u2);

    annotation (Documentation(info="<html>
 
<h2>Provide smooth 2D spline surface and interpolate input data</h2>
This blocks generates a 2D-surface object by approximating spline fitting. 
Data to be fitted may be scattered or on a rectangular grid. 
 
<h3>Rectangular data</h3>
In case of rectangular data points, set boolean parameter rectangular to true and provide the data as array in the common Modelica table form:
 
<pre>
 
  0 | y1 | y2 | y3 | ...
  --+----+----+----+-----
 x1 | z11| z12| z13| ...
 x2 | z21| z22| z23| ...
 x3 | z31| z32| z33| ...
 ...| ...| ...| ...| ...
 
</pre> 
 
<h3>Scattered dat</h3>
In case of scattered data points, set parameter rectangular to false and provide the data as 3- or 4-column array:
 
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
 
 
 
The lasr column w[i] is the optional weight of the data points. If no weights are given, any data points are weighted equally. 
 
</html>"));
  end Surf2d;

end Blocks;
