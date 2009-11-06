within ;
package ApproxSpline "Approximation spline library"
  annotation (uses(Modelica(version="3.0")),
  classOrder = {"Blocks", "Curve1d", "Surf2d", "Plot", "TestModels", "Utilities", "Icons"},
Documentation(info="<html>
<h2> ApproxSpline -- A Modelica library for piecewise polynomial spline approximating</h2>
 
<center>
<img src=\"../images/curve1d.png\" width=350>
<img src=\"../images/surf2d.png\"  width=400>
</center>
 
<p>
This library provides smooth spline approximations of 1 or 2-D data. 
The smoothness and the order of the generated splines can be changed by parameter settings.
ApproxSpline is an interface to some of the functions inplemented in the FITPACK package by Paul Dierckx from netlib.
 
<h3>Installation </h3>
As ApproxSpline is only an interface, you will need the the DIERCKX Fortran routines as well. 
DIERCKX is availabe at netlib: <a href=\"http://netlib.org/dierckx/\"> http://netlib.org/dierckx/ </a>. 
<p>


 
<pre> 

DIERCKX is a package of Fortran subroutines for calculating smoothing splines for various kinds of
data and geometries, with automatic knot selection. This library is also called FITPACK, but is
independent of the FITPACK library by Alan Cline.
 
Reference
    Paul Dierckx, Curve and Surface Fitting with Splines, Oxford University Press, 1993 
Developer
    Paul Dierckx, Department of Computer Science, K.U. Leuven, Celestijnenlaan 200 A, B-3001, Heverlee, Belgium
    Paul.Dierckx@cs.kuleuven.ac.be
</pre>
 
</html>", revisions="<html>
<b>Copyright &copy; 2009, DLR Institute of Vehicle Concepts</b>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<img src=\"../images/dlr_logo.png\"  width=60>

<pre>
Author of DIERCKX library:
    Paul Dierckx
    Department of Computer Science, K.U. Leuven, 
    Celestijnenlaan 200 A, 
    B-3001, Heverlee, 
    Belgium
    email: Paul.Dierckx@cs.kuleuven.ac.be 
</pre>

<dl>
 
  <dt>2009-01-15 <a href='mailto:joerg.ungethuem@dlr.de'> Joerg Ungethuem </a> 
  <dd>realized
 
</dl>
 
</html>"));
end ApproxSpline;
