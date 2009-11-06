within ApproxSpline;
package Utilities "utility functions needed in the library"
  function vmin "return minimum element of vector"
    extends Modelica.Icons.Function;
    input Real vec[:] "input vector";
    output Real y "minimum vector element";
  algorithm
    y :=vec[1];
    for i in 2:size(vec,1) loop
      if vec[i]<y then
        y :=vec[i];
      end if;
    end for;
  end vmin;

  function vmax "return maximum element of vector"
    extends Modelica.Icons.Function;
    input Real vec[:] "input vector";
    output Real y "maximum vector element";
  algorithm
    y :=vec[1];
    for i in 2:size(vec,1) loop
      if vec[i]>y then
        y :=vec[i];
      end if;
    end for;
  end vmax;
end Utilities;
