within ApproxSpline;
package Plot "Utility functions for data visualization"

  model plotDataSurf2d "create data structs to be potted in MATLAB"
  /*
  to plot surface in MATLAB use this command:
 
  load surf2d;figure;hold on;contour3(x,y,z);meshc(x,y,z);xlabel('X');ylabel('Y');hold on;plot3(xd,yd,zd,'*');grid on;view(-50,50);
 
*/

    parameter Boolean rectangular =  false
      "set true if input data is on rectangular grid, otherwise assume scattered data";
    parameter Integer n=20 "number of data points to plot";
    parameter String filename = "surf2d.mat"
      "name of MATLAB formated output file ";
    parameter Real data[ :,:]
      "scattered data to be fitted: [x,y,z,w] (w is optional weight)";
    parameter Real s(min=0) = 0 "smoothing factor (s>0)";

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

    Real[2,1] nk "number of knots in both directions";

  protected
      ApproxSpline.Surf2d.Type spl = ApproxSpline.Surf2d.Type(
        rectangular,
        data,
        s,
        kx,
        ky,
        x_lim,
        y_lim,
        tx,
        ty);

      Integer i;
      Real[n,n] x;
      Real[n,n] y;
      Real[n,n] z;

      Real[if rectangular then (size(data,1)-1)*(size(data,2)-1) else size(data,1),1] xd;
      Real[if rectangular then (size(data,1)-1)*(size(data,2)-1) else size(data,1),1] yd;
      Real[if rectangular then (size(data,1)-1)*(size(data,2)-1) else size(data,1),1] zd;

  algorithm
      for ix in 1:n loop
        for iy in 1:n loop
          x[ix,iy] :=x_lim[1] + (ix - 1)/(n - 1)*(x_lim[2] - x_lim[1]);
          y[ix,iy] :=y_lim[1] + (iy - 1)/(n - 1)*(y_lim[2] - y_lim[1]);
          z[ix,iy] :=ApproxSpline.Surf2d.eval(
            spl,
            x[ix, iy],
            y[ix, iy]);
        end for;
      end for;

      if rectangular then
        i := 0;
        for ix in 2:size(data,1) loop
          for iy in 2:size(data,2) loop
            i :=i + 1;
            xd[i,1] := data[ix, 1];
            yd[i,1] := data[1, iy];
            zd[i,1] := data[ix, iy];
          end for;
        end for;
      else
        i :=0;
        for i in 1:size(data,1) loop
          xd[i,1] :=data[i, 1];
          yd[i,1] :=data[i, 2];
          zd[i,1] :=data[i, 3];
        end for;
      end if;

      nk[1,1] :=1.0*ApproxSpline.Surf2d.getNumberOfKnotsX(spl);
      nk[2,1] :=1.0*ApproxSpline.Surf2d.getNumberOfKnotsY(spl);

      DataFiles.writeMATmatrix(filename, "data", data, false);
      DataFiles.writeMATmatrix(filename, "x", x, true);
      DataFiles.writeMATmatrix(filename, "y", y, true);
      DataFiles.writeMATmatrix(filename, "z", z, true);
      DataFiles.writeMATmatrix(filename, "xd", xd, true);
      DataFiles.writeMATmatrix(filename, "yd", yd, true);
      DataFiles.writeMATmatrix(filename, "zd", zd, true);
      DataFiles.writeMATmatrix(filename, "nk", nk, true);

  end plotDataSurf2d;

  model plotDataCurve1d "create data structs to be potted in MATLAB"
  /*
  to plot curve in MATLAB use this command:
 
  load curve1d;plot(x,y);xlabel('X');ylabel('Y');hold all;plot(xd,yd,'*');hold off;grid on;
 
*/

    parameter Integer n=20 "number of data points to plot";
    parameter String filename = "curve1d.mat"
      "name of MATLAB formated output file ";

    parameter Real data[:,:]
      "data to be fitted: [x,y,w] (w is optional weight)";
    parameter Real s(min=0) = 0 "smoothing factor (s>0)";
    parameter Integer k(min=1,max=5) = 3 "degree of spline (should be odd)";
    parameter Boolean periodic = false
      "set true to generate periodic spline curve";
    parameter Real x_lim[2] = { ApproxSpline.Utilities.vmin(data[:,1]), ApproxSpline.Utilities.vmax(data[:,1])}
      "boundaries of interpolation intervall";
    parameter Real t[:] = fill(0,0)
      "array of knot positions (if given, a least square spline is generated)";

    Real[1,1] nk "number of knots";

  protected
      ApproxSpline.Curve1d.Type spl = ApproxSpline.Curve1d.Type(data, s, k, periodic,x_lim, t);

      Real[n,1] x;
      Real[n,1] y;
      Real[size(data,1),1] xd;
      Real[size(data,1),1] yd;

  algorithm
      for ix in 1:n loop
          x[ix,1] :=x_lim[1] + (ix - 1)/(n - 1)*(x_lim[2] - x_lim[1]);
          y[ix,1] :=ApproxSpline.Curve1d.eval(
            spl,
            x[ix, 1]);
      end for;

      for i in 1:size(data,1) loop
        xd[i,1] :=data[i, 1];
        yd[i,1] :=data[i, 2];
      end for;

      nk[1,1] :=1.0*ApproxSpline.Curve1d.getNumberOfKnots(spl);

      DataFiles.writeMATmatrix(filename, "data", data, false);
      DataFiles.writeMATmatrix(filename, "x", x, true);
      DataFiles.writeMATmatrix(filename, "y", y, true);
      DataFiles.writeMATmatrix(filename, "xd", xd, true);
      DataFiles.writeMATmatrix(filename, "yd", yd, true);
      DataFiles.writeMATmatrix(filename, "nk", nk, true);

  end plotDataCurve1d;

end Plot;
