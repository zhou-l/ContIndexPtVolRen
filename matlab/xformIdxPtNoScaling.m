function Pt = xformIdxPtNoScaling(c1, c2, c3)
if c2 == 0
    u = 0;
    v = -c3/c1;
    if c1 == 0
        v = 0;
    end
    Pt = [u v];
    return;
end

if c1 == -c2
    if rand() < 0.5
        u = -0.5;
    else
        u = 1.5;
    end
    v = 0.0;
    Pt = [u v];
    return;
end
a = -c1/c2;
b = 2 * c3 / (c1 - c2);
theta = atan(a);
if theta > pi/4
    u = theta*2/pi - 1;
else
    u = theta*2/pi + 1;
end
v = (u - 0.5) * b;
Pt = [u v];
end


% void XPCPdataAnalyzer::xform_idxpt_complete_line_desc_noScaling(double c1, double c2, double c3, double& xout, double& yout)
% {
% 
% 	if (c2 == 0) // +/- pi/2; at x1 axis
% 	{
% 		xout = 0.0;
% 		yout = -c3 / c1;
% 		if (c1 == 0.0)
% 			yout = 0.0;
% 		return;
% 	}
% 
% 	if (c1 == -c2) // pi/4
% 	{
% 		xout = (g_randGen.rand() < 0.5) ? -0.5 : 1.5; // half chance to the left, half to the right
% 		yout = 0.0;
% 		return;
% 	}
% 
% 	// parameters of the associated Cartesian line 
% 	double a = -c1 / c2;
% 	double b = 2 * c3 / (c1 - c2);  //-c3 / c2;
% 
% 	double theta = atan(a);
% 	double scale = 1.0;
% 
% 	if (theta > M_PI / 4) // left of the axes pair
% 	{
% 		xout = theta * 2 / M_PI - 1;
% 	}
% 	else // For theta < 0 && theta < pi/4
% 	{
% 		xout = theta * 2 / M_PI + 1;
% 	}
% 	// Old transformation method with the origin at 0
% 	//yout = xout * b * scale; // output y
% 	// with new offset
% 	yout = (xout - 0.5f) * b * scale; // output y
% }