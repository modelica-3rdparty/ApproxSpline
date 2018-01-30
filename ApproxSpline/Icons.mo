within ApproxSpline;
package Icons
  extends Modelica.Icons.IconsPackage;
  class External "Icon of an external object"
      annotation(Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{100,100}}), graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            fillColor={181,181,181},
            fillPattern=FillPattern.Solid), Text(
            extent={{-94,94},{94,-94}},
            fillColor={181,181,181},
            fillPattern=FillPattern.Solid,
            textString="X")}),
                        Documentation(info="<html>
<p>
This icon is designed for an <b>ExternalObject</b> type.
</p>
</html>"));

  end External;

  class Curve1d

    annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{100,100}}), graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            fillColor={250,250,250},
            fillPattern=FillPattern.Solid), Bitmap(extent={{-96,96},{96,-96}},
              fileName="modelica://ApproxSpline/Resources/Images/curve1d.png")}));
  end Curve1d;

  class Surf2d

    annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}), graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            fillColor={250,250,250},
            fillPattern=FillPattern.Solid), Bitmap(extent={{-98,98},{98,-98}},
              fileName="modelica://ApproxSpline/Resources/Images/surf2d.png")}));
  end Surf2d;
end Icons;
