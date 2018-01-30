within ;
package ApproxSpline "Approximation spline library"
  extends Modelica.Icons.Package;
  annotation (uses(Modelica(version="3.2.2")), version="1.0.0",
Documentation(info="<html>
<h4>ApproxSpline -- A Modelica library for piecewise polynomial spline approximating</h4>

<center>
<img src=\"modelica://ApproxSpline/Resources/Images/curve1d.png\" width=350>
<img src=\"modelica://ApproxSpline/Resources/Images/surf2d.png\"  width=400>
</center>

<p>
This library provides smooth spline approximations of 1 or 2-D data.
The smoothness and the order of the generated splines can be changed by parameter settings.
ApproxSpline is an interface to some of the functions inplemented in the FITPACK package by Paul Dierckx from netlib.

<h5>Installation</h5>
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
<img src=\"modelica://ApproxSpline/Resources/Images/dlr_logo.png\"  width=60>

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

  <dt>2009-01-15 <a href=\"http://www.dlr.de/tt/desktopdefault.aspx/tabid-4074/6449_read-26380\">J&ouml;rg Weiss-Ungeth&uuml;m</a>
  <dd>realized

</dl>

</html>"));
end ApproxSpline;
