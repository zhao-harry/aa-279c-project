function Mc = controlLaw(alpha,alphadot,control)
    Mc = -control.Kp .* alpha - control.Kd .* alphadot;
end