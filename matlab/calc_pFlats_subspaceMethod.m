function idxPt = calc_pFlats_subspaceMethod(mu, majEig, secEig)
N = size(mu,1);
ptMajEig = majEig+mu;

end

% int N = mu.size(); // dimension of the input data
% 
% 	Eigen::VectorXf ptMajEig = majEig + mu;
% 	vector<float> vPtMajEig(ptMajEig.data(), ptMajEig.data() + ptMajEig.rows());
% 	vector<float> vPtMu(mu.data(), mu.data() + mu.rows());
% 	// 0. Set N-D value
% 	xpcp_tuple.x = vPtMu;
% 	// threshold the length of the major eigen vector
% 	float twoFlat_discard_thres = 1e-3f;
% 	float oneFlat_discard_thres = 1e-3f;
% 
% 	// 1. Compute 1-flat
% 	vector<Point> pcp_v0_no_repeat;
% 
% 	if (g_params.XformMethod() == XF_NGUYEN_ROSEN)
% 	{
% 		_pcp->calcOneFlatParametricForm(vPtMu, vPtMajEig, pcp_v0_no_repeat);
% 	}
% 	else
% 		_pcp->calcOneFlatGeneralForm(vPtMu, vPtMajEig, pcp_v0_no_repeat);
% 	xpcp_tuple.pcp_1flat.resize(pcp_v0_no_repeat.size());
% 
% 	double u = 0.0;
% 	double v = 0.0;
% 	// for nguyen & rosen's tvcg 2017
% 	////////////////////////////////////
% 	double uu = 0.0;
% 	double vv = 0.0;
% 	double x0 = 0.0;
% 	double y0 = 0.0;
% 	bool isRotated = false;
% 	///////////////////////////////////
% 	for (size_t i = 0; i < pcp_v0_no_repeat.size(); i++)
% 	{
% 		Point idxPt = pcp_v0_no_repeat[i]; // 3D line coordinates0
% 
% 		if (g_params.XformMethod() == XF_IDX_ONLY || g_params.XformMethod() == XF_IDX_LINE)
% 		{
% 			double c1 = idxPt.x;
% 			double c2 = idxPt.y;
% 			double c3 = idxPt.z;
% 			if (idxPt.y != 0.0)
% 			{
% 				double m = -idxPt.x / idxPt.y;
% 				double b = -idxPt.z / idxPt.y;
% 
% 				c1 = -m;
% 				c2 = 1.0;
% 				c3 = -b;
% 			}
% 
% 			if(g_params.IsScaleXform())
% 				xform_idxpt_complete_line_desc(c1, c2, c3, u, v);
% 			else
% 				xform_idxpt_complete_line_desc_noScaling(c1, c2, c3, u, v);
% 		}
% 		else if (g_params.XformMethod() == XF_NGUYEN_ROSEN)
% 		{
% 			x0 = vPtMu[i];
% 			y0 = vPtMu[i + 1];
% 			uu = idxPt.x;
% 			vv = idxPt.y;
% 			xform_nguyen_rosen_tvcg(uu, vv, x0, y0, u, v, isRotated);
% 		}
% 		else
% 			lineCoord2pointCoord2D(idxPt.x, idxPt.y, idxPt.z, u, v);
% 		// Adding the offset of the dimension
% 		xpcp_tuple.pcp_1flat[i] = FLOATVECTOR2(float(u) + float(i), float(v));
% 		float eigVLen = Eigen::Vector2f(majEig[i], majEig[i + 1]).norm();
% 		///////////////////////////////////////////////////////////
% 		// Record strength of the eigen vector in each subspace and get the min/max
% 		xpcp_tuple.strength_1flat[i].x = eigVLen;
% 		if (g_params.XformMethod() == XF_NGUYEN_ROSEN) //HACK: use the y component of the strength to encode whether the point is rotated for Nguyen&Rosen!!!
% 			xpcp_tuple.isRotated_1flat[i] = isRotated;
% 		//**********************************************************
% 		// Alternatively, use the correlation value for 1-flat strength
% 		// NOTE: Uncomment to use correlation value as strength!
% 		//	xpcp_tuple.strength_1flat[i].x = abs(corr_mat(i,i+1));
% 		//**********************************************************
% 		float minEigVlen = MIN(_1flatsMinMaxPerSubspace[i].x, xpcp_tuple.strength_1flat[i].x);
% 		_1flatsMinMaxPerSubspace[i].x = minEigVlen;
% 		float maxEigVlen = MAX(_1flatsMinMaxPerSubspace[i].y, xpcp_tuple.strength_1flat[i].x);
% 		_1flatsMinMaxPerSubspace[i].y = maxEigVlen;
% 	}
% 
% 	//=======================================================================			
% 	// 2. Compute 2-flat
% 	// For now, use 3-tuple in sequency of the axis ordering
% 	vector<FLOATVECTOR2> indexPt;
% 	bool oneIndex = !g_use_repeat_pcp;
% 	int numIdxPerPlane = oneIndex ? 1 : 4;
% 	xpcp_tuple.pcp_2flat.resize((N - 2) * numIdxPerPlane);
% 	int normMethod = 1; // 0: length of the normal vector, 1: product of the lengths of two tangent vectors, 2: the larger of the tangents, 3: the smaller of the the tangents
% 	for (int j = 0; j < N - 2; j++)
% 	{
% 		// Get 3D subspace!!! Use projected eigenvectors as tangent vectors of the plane
% 		Eigen::Vector3f mu3d = mu.block(j, 0, 3, 1);
% 		Eigen::Vector3f majEig3d = majEig.block(j, 0, 3, 1);
% 		Eigen::Vector3f secEig3d = secEig.block(j, 0, 3, 1);
% 		//majEig3d.normalize(); // normalize these tangent vectors
% 		//secEig3d.normalize();
% 
% 		// Do analytical computation to find 2-flat for the subspace
% 		Eigen::Vector3f norm3d = majEig3d.cross(secEig3d);
% 		float strength_2flat = 0.0f;
% 		if (norm3d == Eigen::Vector3f(0, 0, 0))
% 		{
% 			xpcp_tuple.strength_2flat[j].x = strength_2flat;
% 			xpcp_tuple.strength_2flat[j].y = strength_2flat;
% 
% 			_2flatsMinMaxPerSubspace[j].x = MIN(strength_2flat, _2flatsMinMaxPerSubspace[j].x);
% 			_2flatsMinMaxPerSubspace[j].y = MAX(strength_2flat, _2flatsMinMaxPerSubspace[j].y);
% 			continue;
% 		}
% 
% 		/////////////////////////////////////
% 
% 		// Normalize the normal vector
% 		strength_2flat = computeNormalVecStrength(majEig3d, secEig3d, norm3d, normMethod);
% 		//strength_2flat = norm3d.norm();
% 		//////////////////////////////////////////
% 
% 		_2flatsMinMaxPerSubspace[j].x = MIN(strength_2flat, _2flatsMinMaxPerSubspace[j].x);
% 		_2flatsMinMaxPerSubspace[j].y = MAX(strength_2flat, _2flatsMinMaxPerSubspace[j].y);
% 
% 		PlaneCoeffs Pl;
% 		// The plane needs a normalized normal vector
% 		norm3d.normalize();
% 		Pl.c1 = norm3d.x();
% 		Pl.c2 = norm3d.y();
% 		Pl.c3 = norm3d.z();
% 		Pl.c0 = norm3d.x() * mu3d[0] + norm3d.y() * mu3d[1] + norm3d.z() * mu3d[2];
% 
% 		// Compute indexed point(s)
% 		_pcp->calc2FlatFrom3DPlane(Pl, indexPt, oneIndex);
% 
% 
% 		// Set 2-flat
% 		for (int k = 0; k < numIdxPerPlane; k++)
% 		{
% 			// Conduct the same transformation as 1-flats
% 
% 			if (g_params.XformMethod() == XF_IDX_ONLY || g_params.XformMethod() == XF_IDX_LINE)
% 			{
% 				double u = 0.0;
% 				double v = 0.0;
% 				xform(indexPt[k].x, indexPt[k].y, j, u, v);
% 				indexPt[k].y = v;
% 				indexPt[k].x = u; // offset by 0.5(?) to match the width of axes.
% 			}
% 
% 			// Shift to the starting dimension
% 
% 
% 			indexPt[k].x += float(j) + 0.5f;
% 			xpcp_tuple.pcp_2flat[numIdxPerPlane * j + k] = indexPt[k];
% 			xpcp_tuple.strength_2flat[j].x = strength_2flat; // set strength
% 		}
% 
% 	}

% void PCPInselberg::calcOneFlatGeneralForm(const vector<float>& cart_pm1flat_1, const vector<float>& cart_pm1flat_2, std::vector<Point>& pcp_oneflat)
% {
% 	pcp_oneflat.resize(cart_pm1flat_1.size() - 1);
% 	for (size_t i = 0; i < cart_pm1flat_1.size() - 1; i++)
% 	{
% 		// compute the intersections formed by line segements connecting neighboring points of p-1 flat 1 and p-1 flat 2
% 
% 		// First 2D cartesian point 
% 		Point p1(cart_pm1flat_1[i], cart_pm1flat_1[i+1]);
% 	
% 		// Second p-1 flat
% 		Point p2(cart_pm1flat_2[i], cart_pm1flat_2[i + 1]);
% 
% 		// general form of the line
% 		float c1 = p2.y - p1.y;
% 		float c2 = -(p2.x - p1.x);
% 		float c3 = p1.y * (p2.x - p1.x) - p1.x * (p2.y - p1.y);
% 		
% 		pcp_oneflat[i] = Point(c1, c2, c3);
% 	}
