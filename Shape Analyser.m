(* ::Package:: *)

BeginPackage["ShapeanAlyser`"]


(*Utils*)
launchNukes::usage = "launchNukes[] launches 2*$ProcessorCount-1 - Length@Kernels[] Kernels"


(*Get Phase Boundary Points & Bounds*)
getPixelCoordinates::usage = "getPixelCoordinates[image] finds boundary of the mask"
getBoundingBox::usage = "getBoundingBox[image] finds Bounding Box of the mask"
constructBubbleRegion::usage = "constructBubbleRegion[boundaryPoints] constructs MeshRegion encompassed by boundaryPoints"


(*Curve evolution*)
kc1::usage= "kc1[a,b,c,d,e,f]"
kc2::usage= "kc2[a,b,c,d,e,f]"
curveEvolution::usage= "curveEvolution[polygon, fraq, relevanceFunction]"


(*Curvature analasys*)
getCurvatureProfile::usage = "getCurvatureProfile[boundaryPoints, windowFactor]"
getCurvatureExtrema::usage = "getCurvatureExtrema[curvatureProfile]"
mapExtremaToBoundary::usage = "mapExtremaToBoundary[boundarySpline, curvatureProfile, curvatureExtrema] ... if no curvatureExtreama specified, all will be used"


(*Boundary refinements*)
buildBoundarySpline::usage = "buildBoundarySpline[boundaryPoints, splineOrder] builds cosed spline of order splineOrder for boundaryPoints"
uniformlySampleSpline::usage = "uniformlySampleSpline[spline, samplingPoints] generates a number of samplingPoints of points uniformly spaced along the spline" 
downsampleBoundary::usage = "downsampleBoundary[boundaryPoints, boundarySpline, tolerance, windowFactor]"


(*Voronoi Shiz*)
getVoronoiInteriorNodes::usage = "getVoronoiInteriorNodes[boundaryPoints, bubbleRegion]"
getAssistingSkeleton::usage = "getAssistingSkeleton[testImage]"
getVoronoiSpanningTree::usage = "getVoronoiSpanningTree[voronoiInteriorNodes]"
getVoronoiDendrites::usage =  "getVoronoiDendrites[voronoiSpanningTree]"


(*Camber line construction*)
constructCamberLine::usage =  "constructCamberLine[voronoiInteriorNodes, voronoiSpanningTree]"
finishCamberLine::usage = "finishCamberLine[camberMain, extremaLocations]"
resampleCamberLine::usage = "resampleCamberLine[camberLine, splineOrder, samplingPoints]"
getChord::usage = "getChord[camberLine]"


(*Measurements*)
getChordParams::usage = "getChordParams[chord]"
getCamberLineParams::usage = "getCamberLineParams[camberLine]"


Begin["Private`"]


useKernels=2*$ProcessorCount-1;

launchNukes[]:=If[
Length@Kernels[]<useKernels,
Speak@"Launching Nukes";
LaunchKernels[
useKernels-Length@Kernels[]
],
Nothing
];


getPixelCoordinates[image_]:=DeleteDuplicates@First@List@Delete[#,0]&@First@Last@First@ComponentMeasurements[image,"Contours"];


getBoundingBox[image_]:=Transpose@First@{1}/.ComponentMeasurements[image,"BoundingBox"];


constructBubbleRegion[boundaryPoints_]:=Module[{boundaryLines},

boundaryLines=AppendTo[
Table[
Line@{k,k+1},
{k,1,Length@boundaryPoints-1}
],
Line@{Length@boundaryPoints,1}
];//Quiet;

BoundaryMeshRegion[boundaryPoints,boundaryLines]
]


buildBoundarySpline[boundaryPoints_, splineOrder_]:=BSplineFunction[
boundaryPoints,
SplineDegree->splineOrder,SplineClosed->True
];


arcIncrement[t_?NumericQ,spline_]:=If[
t-1.<=0,
spline'[t],
spline'[1]
]


getCheckpoints[spline_]:=NDSolveValue[

{
t'[s]==1/Norm@arcIncrement[t[s],spline],
t[0]==0,
WhenEvent[t[s]==1,"StopIntegration"]
},

t,{s,0,1+NIntegrate[
Norm[spline'[t]],
{t,0,1},Method->"LocalAdaptive",MaxRecursion->500
]}

];


uniformlySampleSpline[spline_,pointCount_]:=Module[{arcCheckpoints},
arcCheckpoints=getCheckpoints@spline;
Map[
spline,

arcCheckpoints@Rescale[
Range[0,1,1/pointCount],
{0,1},
First@arcCheckpoints["Domain"]
]

]
];


interpolateAndUniformlyUpsampleBoundary[boundaryPoints_, splineOrder_, samplingPoints_]:=Module[{spline},
spline = buildBoundarySpline[boundaryPoints, splineOrder];
uniformlySampleSpline[spline,samplingPoints]
];


\[Beta]c=Compile[{{a,_Real},{b,_Real},{c,_Real},{d,_Real},{e,_Real},{f,_Real}},With[{e1={c,d}-{a,b},e2={e,f}-{c,d}},ArcCos[e1.e2/(Norm[e1] Norm[e2])]]];

kc1=Compile[{{a,_Real},{b,_Real},{c,_Real},{d,_Real },{e,_Real},{f,_Real}},With[{n1=Norm[{c,d}-{a,b}],n2=Norm[{e,f}-{c,d}]},(\[Beta]c[a,b,c,d,e,f]n1 n2)/(n1+n2)]];

kc2=Compile[{{a,_Real},{b,_Real},{c,_Real},{d,_Real },{e,_Real},{f,_Real}},With[{n1=Norm[{c,d}-{a,b}],n2=Norm[{e,f}-{c,d}]},(Abs[\[Beta]c[a,b,c,d,e,f]-\[Pi]]n1 n2)/(n1+n2)]];

reductionStep[polygon_, relevanceFunction_]:=Module[{segments, scores, worstPoint},
segments=Line@Range@Length@polygon;
scores=Map[relevanceFunction[Sequence@@(Flatten@Part[polygon,#])]&,Partition[First@segments,3,1,{2,2}]]//Quiet;
worstPoint=First@Ordering[scores];
Delete[polygon,worstPoint]
];

curveEvolution[polygon_,fraq_,relevanceFunction_]:=Nest[reductionStep[#,relevanceFunction] &,polygon,Round[fraq*Length[polygon]]];


getCurvatureProfile[boundaryPoints_,windowFactor_]:=Module[{adaptiveWindow,dx1,dy1,dx2,dy2},
adaptiveWindow=If[#<1,1,#]&@(windowFactor*Length@boundaryPoints);
{{dx1,dy1},{dx2,dy2}}=Table[
GaussianFilter[
boundaryPoints[[All,dimensions]],
adaptiveWindow,derivatives
],
{derivatives,2},
{dimensions,2}
];
(dx1*dy2-dx2*dy1)/((dx1^2+dy1^2)^(3/2))
];


getCurvatureExtrema[curvatureProfile_]:=Select[
FindPeaks@curvatureProfile,
Last@#>0&
];


mapExtremaToBoundary[spline_, curvatureProfile_]:=Module[{curvatureExtrema, splineParameters, criticalParameters},

curvatureExtrema=getCurvatureExtrema[curvatureProfile];

splineParameters=N@Rescale@Range@Length@curvatureProfile;
criticalParameters=splineParameters[[
curvatureExtrema[[;;,1]]
]];
spline/@criticalParameters
]

mapExtremaToBoundary[spline_, curvatureProfile_, curvatureExtrema_]:=Module[{splineParameters, criticalParameters},

splineParameters=N@Rescale@Range@Length@curvatureProfile;
criticalParameters=splineParameters[[
curvatureExtrema[[;;,1]]
]];
spline/@criticalParameters
]


dist[q:{x_,y_}, {p1:{x1_,y1_},p2:{x2_,y2_}}] := With[

{u =(q-p1).(p2-p1)/(p2-p1).(p2-p1)},

Which[
u<=0, Norm[q-p1],
u>=1, Norm[q-p2],
True, Norm[q-(p1+u(p2-p1))]
]

];

testSeg[seg[points_List], tol_] := Module[

{},

dists = dist[#,{
points[[1]],points[[-1]]
}]& /@ points[[
Range[2,Length[points]-1]
]];

max = Max@dists;

If[
max > tol,
pos = Position[dists,max][[1,1]]+1;
{
seg@points[[Range[1,pos]]],
seg@points[[
Range[pos,Length@points]
]]
},
seg[points,done]
]
] /; Length@points> 2;

testSeg[seg[points_List], tol_] := seg[points, done];

testSeg[seg[points_List,done], tol_] := seg[points,done];

dpSimp[points_, tol_] := Append[First /@ First /@ Flatten[{seg@points} //. s_seg :> testSeg[s, tol]],Last[points]]


downsampleBoundary[upsampledBoundary0_, spline_, tolerance_, windowFactor_]:=Module[
{upsampledBoundary,curvatureExtrema,curvatureProfile,
topTwoExtrema,topTwoExtremaLocations,
pointA, pointB,
positionA, positionB,
segmentA, segmentB,
pointsA,pointsB
},

upsampledBoundary=upsampledBoundary0//DeleteDuplicates;

curvatureProfile=getCurvatureProfile[upsampledBoundary, windowFactor];
curvatureExtrema=getCurvatureExtrema[curvatureProfile];

topTwoExtrema = #[[1;;2]]&@ReverseSortBy[#,Last]&@curvatureExtrema;
topTwoExtremaLocations=mapExtremaToBoundary[spline, curvatureProfile, topTwoExtrema];

{pointA, pointB} = First@Nearest[upsampledBoundary, #@topTwoExtremaLocations,1]&/@{First,Last};
{positionA,positionB} = First@First@Position[upsampledBoundary,#]&/@{pointA, pointB};

If[
positionA >positionB,


segmentA=Join[
upsampledBoundary[[positionA;;]],
upsampledBoundary[[1;;positionB]]
];

segmentB=upsampledBoundary[[positionB;;positionA]];,


segmentA=upsampledBoundary[[positionA;;positionB]];

segmentB=Join[
upsampledBoundary[[positionB;;]],
upsampledBoundary[[1;;positionA]]
];

];
{pointsA,pointsB}=dpSimp[#,tolerance]&/@{segmentA,segmentB};
Join[pointsA,pointsB]//DeleteDuplicates
];


getVoronoiInteriorNodes[boundaryPonints_,bubbleRegion_]:=Module[{voronoiPoints},
voronoiPoints=VoronoiMesh@boundaryPonints;
Select[
MeshCoordinates@voronoiPoints,
RegionMember[bubbleRegion,#]==True&
]
];


getAssistingSkeleton[image_]:=Module[{camberSkeleton, skeletonPoints},
camberSkeleton=Thinning@image;
skeletonPoints=getPixelCoordinates@camberSkeleton;
N@Map[
Mean,
GatherBy[skeletonPoints,First]
]
];


getVoronoiSpanningTree[voronoiNodes_]:=Module[{distanceMatrix, getGraph},
distanceMatrix=DistanceMatrix[
voronoiNodes,
DistanceFunction->EuclideanDistance
];
getGraph=WeightedAdjacencyGraph[
distanceMatrix,DirectedEdges->False
];
FindSpanningTree@getGraph
];


getVoronoiDendrites[graphSpanningTree_]:=Position[
VertexDegree@graphSpanningTree,1
][[All,1]];


constructCamberLine[voronoiNodes_, graphSpanningTree_]:=Module[{dendriteEndpoints, allPaths, longestPath},
dendriteEndpoints=getVoronoiDendrites[graphSpanningTree];
allPaths=Flatten[
Parallelize@Outer[
FindShortestPath[graphSpanningTree,#1,#2]&,
dendriteEndpoints,dendriteEndpoints
],1
];

longestPath=First@MaximalBy[allPaths,Length];
voronoiNodes[[longestPath]]
]


finishCamberLine[camberMain_, extremaLocations_]:=Module[{camberEndpoints},
camberEndpoints=Map[
Nearest[extremaLocations,#@camberMain]&,
{First,Last}
];
Flatten[
{
First@camberEndpoints,
camberMain,
Last@camberEndpoints
},1
]
]


resampleCamberLine[camberLine_, splineOrder_,samplingPoints_]:=Module[{camberSpline,arcCheckpoints},
camberSpline=BSplineFunction[
camberLine,
SplineDegree->splineOrder,SplineClosed->False
];

arcCheckpoints=getCheckpoints@camberSpline;

uniformlySampleSpline[camberSpline,samplingPoints]
]


getCamberLineParams[camberLine_]:={"lenght"->
Total@MapThread[

EuclideanDistance[
camberLine[[#1]],camberLine[[#2]]
]&,

{
Range[1,Length@camberLine-1],
Range[2,Length@camberLine]
}

]
}


getChord[camberLine_]:=Map[#@camberLine&,{First,Last}];


getChordParams[chord_]:=Module[{chordLenght,principalChord,tiltAngle},
chordLenght=EuclideanDistance[First@chord, Last@chord];
principalChord[x_]=Normal@LinearModelFit[chord,x,x];
tiltAngle=(180/\[Pi])*ArcTan[principalChord'[x]];
{"lenght"->chordLenght, "tilt"->tiltAngle}
]


End[]

EndPackage[]
