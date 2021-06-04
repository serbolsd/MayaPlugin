//-
// ==========================================================================
// Copyright 1995,2006,2008 Autodesk, Inc. All rights reserved.
//
// Use of this software is subject to the terms of the Autodesk
// license agreement provided at the time of installation or download,
// or which otherwise accompanies this software in either electronic
// or hard copy form.
// ==========================================================================
//+

// Example Plugin: circleNode.cpp
//
// This plug-in is an example of a user-defined dependency graph node.
// It takes a number as input (such as time) and generates two output
// numbers one which describes a sine curve as the input varies and
// one that generates a cosine curve. If these two are hooked up to
// the x and z translation attributes of an object the object will describe
// move through a circle in the xz plane as time is changed.
//
// Executing the command "source circleNode" will run the MEL script which will
// create a new "Circle" menu with a single item. Selecting this will build
// a simple model (a sphere which follows a circular path) which can be played back,
// by clicking on the "play" icon on the time slider.  Note: the circleNode
// plugin needs to be loaded before the "Circle" menu item can be executed
// properly.
//
// The node has two additional attributes which can be changed to affect
// the animation, "scale" which defines the size of the circular path, and
// "frames" which defines the number of frames required for a complete circuit
// of the path. Either of these can be hooked up to other nodes, or can
// be simply set via the MEL command "setAttr" operating on the circle node
// "circleNode1" created by the MEL script. For example, "setAttr circleNode1.scale #"
// will change the size of the circle and "setAttr circleNode1.frames #"
// will cause objects to complete a circle in indicated number of frames.



#include <maya/MPxNode.h> 

#include <maya/MString.h> 
#include <maya/MTypeId.h> 
#include <maya/MPlug.h>

#include <maya/MFnNumericAttribute.h>
#include <maya/MVector.h>

#include <maya/MFnPlugin.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>

#include <string.h> 
#include <sys/types.h>
#include <maya/MStatus.h>
#include <maya/MPxCommand.h>
#include <maya/MStringArray.h>
#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MFnMesh.h>
#include <maya/MFnSet.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItMeshEdge.h>
#include <maya/MFloatVector.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFloatArray.h>
#include <maya/MObjectArray.h>
#include <maya/MObject.h>
#include <maya/MPxFileTranslator.h>
#include <maya/MFnDagNode.h>
#include <maya/MItDag.h>
#include <maya/MDistance.h>
#include <maya/MIntArray.h>
#include <maya/MIOStream.h>
#include <maya/MDrawProcedureBase.h>

#include <math.h>

#include <maya/MSimple.h>

#include <vector>

namespace TYPESAVED {
		enum E {
				MODEL = 0,
				TEXTURE
		};
}

struct FILE_BASE_INFO {
		unsigned int id;
		unsigned int NameSize;
		unsigned int DirectionSize;
		TYPESAVED::E type;
};

struct NAME_INFO {
		std::string name;
		std::string dir;
};

struct VertexData {
		MPoint vertex;
		MVector normal;
		MVector tanget;
		MVector binormal;
		float2 uv;
};

class ExportComand : public MPxCommand {
public:
		ExportComand() {};
		virtual MStatus	doIt(const MArgList&);
		static void* creator();
		void extractVertices();
		void extractIndex();
 
		std::vector<VertexData> m_mesh;
};																

MStatus ExportComand::doIt(const MArgList&)
{
		extractIndex();
		//extractVertices();
		return MS::kSuccess;
}

void* ExportComand::creator()
{																
	return new ExportComand;
}

void ExportComand::extractVertices()
{
		// we assume here that Maya has been initialized and the file in
		// question has already been loaded.

		MStatus stat;
		MItDag dagIter(MItDag::kBreadthFirst, MFn::kInvalid, &stat);

		for (; !dagIter.isDone(); dagIter.next())
		{
				MDagPath dagPath;
				stat = dagIter.getPath(dagPath);

				if (stat)
				{
						MFnDagNode dagNode(dagPath, &stat);

						// this object cannot be intermediate, and it must be a mesh
						// and it can't be a transform.
						// Intermediate objects are special meshes
						// that are not drawn used for mesh morphs or something.
						if (dagNode.isIntermediateObject()) continue;
						if (!dagPath.hasFn(MFn::kMesh)) continue;
						if (dagPath.hasFn(MFn::kTransform)) continue;

						MFnMesh fnMesh(dagPath);

						// get the vertices that are part of the current mesh
						MPointArray vertexList;
						MFloatVectorArray normalList;
						MFloatVectorArray tangentList;
						MFloatVectorArray biNorList;
						MFloatArray uArray, vArray;

						fnMesh.getPoints(vertexList, MSpace::kWorld);
						fnMesh.getNormals(normalList, MSpace::kWorld);
						fnMesh.getTangents(tangentList, MSpace::kWorld);
						fnMesh.getBinormals(biNorList, MSpace::kWorld);
						fnMesh.getUVs(uArray, vArray);
						//fnMesh.gettangen
						//fnMesh.getnormal
						//fnMesh.tanget

						int vertexLength = vertexList.length();
						int normalLength = normalList.length();
						int tangentLength = tangentList.length();
						int biNorLength = biNorList.length();
						int ULength = uArray.length();
						int VLength = vArray.length();
						// iterate through all the vertices

						for (INT32 i = 0; i < vertexList.length(); i++)
						{
								vertexList[i].cartesianize();
								MPoint point = vertexList[i];

								// here is your data... now go do whatever you want with
								// it. if you need a unique identifier for this vertex,
								// use it's index in the mesh, and some kind of mesh id.
								// these stay constant while exporting ( so long as the file is
								// not edited )
								(point.x, point.y, point.z);
						}
				}
		}
}

void ExportComand::extractIndex() {
  MStatus stat;
  MItDag dagIter(MItDag::kBreadthFirst, MFn::kInvalid, &stat);

  for (; !dagIter.isDone(); dagIter.next())
  {
    MDagPath dagPath;
    stat = dagIter.getPath(dagPath);

    if (stat)
    {
      MFnDagNode dagNode(dagPath, &stat);

      if (dagNode.isIntermediateObject()) continue;
      if (!dagPath.hasFn(MFn::kMesh)) continue;
      if (dagPath.hasFn(MFn::kTransform)) continue;

      // get the mesh and all its vertices
      MFnMesh fnMesh(dagPath);
      MPointArray vertexList;
      fnMesh.getPoints(vertexList, MSpace::kObject);
						MFloatVectorArray tangentList;
						fnMesh.getTangents(tangentList, MSpace::kObject);
						MFloatVectorArray normalList;
						fnMesh.getNormals(normalList, MSpace::kObject);

						int vertexLength = vertexList.length();
						int normalLength = normalList.length();
						int tangentLength = tangentList.length();
						std::vector<bool> vertexPerIndex(false);
						vertexPerIndex.resize(vertexList.length());
						m_mesh.clear();
      MObject comp;
      // now iterate over all the polygons in the mesh
      MItMeshPolygon piter(dagPath, comp);
      for (; !piter.isDone(); piter.next())
      {
								
								int numTriangles;
								piter.numTriangles(numTriangles);
								for (size_t i = 0; i < numTriangles; i++)
								{
										MPointArray pA;
										MIntArray iA;
										stat = piter.getTriangle(i, pA, iA); 
										if (!stat)
										{
												continue;
										}
										int indexNum = iA.length();
										int verNum = pA.length();
										for (size_t x = 0; x < indexNum; x++)
										{
												int index = iA[x];
										}
										for (size_t x = 0; x < verNum; x++)
										{
												MPoint index = pA[x];
										}
								}
								MPointArray pA;
								MIntArray iA;
								stat = piter.getTriangles(pA, iA);
								if (!stat)
								{
										continue;
								}
								int numIndex = iA.length();
								for (size_t i = 0; i < iA.length(); i++)
								{
										int index = iA[i];
								}
        // for each polygon you can get the indices of
        // each of its vertices in the vertex list above
        MIntArray vertexIdxs;
        piter.getVertices(vertexIdxs);
								for (int i = 0; i < vertexIdxs.length(); i++)
								{
										int paver = vertexIdxs[i];
										//m_mesh[vertexIdxs[i]].vertex = vertexList[vertexIdxs[i]];
										MPoint vertex = vertexList[vertexIdxs[i]];
										MVector N;
										fnMesh.getFaceVertexNormal(piter.index(), vertexIdxs[i], N);
										//m_mesh[vertexIdxs[i]].normal = N;

										MVector T;
										fnMesh.getFaceVertexTangent(piter.index(), vertexIdxs[i], T);
										//m_mesh[vertexIdxs[i]].tanget = T;

										MVector BN;
										fnMesh.getFaceVertexBinormal(piter.index(), vertexIdxs[i], BN);
										//m_mesh[vertexIdxs[i]].binormal = BN;

										//piter.getUV(vertexIdxs[i], m_mesh[vertexIdxs[i]].uv);
								}
        if (vertexIdxs.length() == 3)
        {
          // process a triangle
          MPoint point0 = vertexList[vertexIdxs[0]];
          MPoint point1 = vertexList[vertexIdxs[1]];
          MPoint point2 = vertexList[vertexIdxs[2]];

          //processTriangle(point0, point1, point2);
        }
      }

						//for (int i = 0; i < vertexLength; i++)
						//{
						//		m_mesh[i].normal.normalize();
						//		m_mesh[i].tanget.normalize();
						//		m_mesh[i].binormal.normalize();
						//}
    }
  }
}


MStatus initializePlugin( MObject _obj )						
{																
	MFnPlugin	plugin( _obj, "Autodesk",  "2020" );
	MStatus		stat;											
	stat = plugin.registerCommand( "ExportComand",					
	                                ExportComand::creator );	    
	if ( !stat )												
		stat.perror("registerCommand");							
	return stat;												
}									

MStatus uninitializePlugin( MObject _obj )						
{																
	MFnPlugin	plugin( _obj );									
	MStatus		stat;											
	stat = plugin.deregisterCommand( "ExportComand" );				
	if ( !stat )												
		stat.perror("deregisterCommand");						
	return stat;												
}

//
//class JDExport {
//public:
//		JDExport() {};
//		virtual         ~JDExport() {};
//		static void* creator();
//		void extractVertices();
//};
//
//// The circle class defines the attributes
//// and methods necessary for the circleNode plugin
////
//class circle : public MPxNode
//{
//public:
//		circle();
//		~circle() override;
//
//		MStatus		compute(const MPlug& plug, MDataBlock& data) override;
//
//		static	void* creator();
//		static	MStatus		initialize();
//
//public:
//		static	MObject		input;		// The input value.
//		static	MObject		sOutput;	// The sinusoidal output value.
//		static	MObject		cOutput;	// The cosinusoidal output value.
//		static	MObject		frames;		// Number of frames for one circle.
//		static	MObject		scale;		// Size of circle.
//		static	MTypeId		id;
//};
//
//MTypeId     circle::id(0x80005);
//MObject     circle::input;
//MObject     circle::sOutput;
//MObject     circle::cOutput;
//MObject	    circle::frames;
//MObject 	circle::scale;
//
//
//void JDExport::extractVertices() {
//		// we assume here that Maya has been initialized and the file in
//		// question has already been loaded.
//
//		MStatus stat;
//		MItDag dagIter(MItDag::kBreadthFirst, MFn::kInvalid, &stat);
//
//		for (; !dagIter.isDone(); dagIter.next())
//		{
//				MDagPath dagPath;
//				stat = dagIter.getPath(dagPath);
//
//				if (stat)
//				{
//						MFnDagNode dagNode(dagPath, &stat);
//
//						// this object cannot be intermediate, and it must be a mesh
//						// and it can't be a transform.
//						// Intermediate objects are special meshes
//						// that are not drawn used for mesh morphs or something.
//						if (dagNode.isIntermediateObject()) continue;
//						if (!dagPath.hasFn(MFn::kMesh)) continue;
//						if (dagPath.hasFn(MFn::kTransform)) continue;
//
//						MFnMesh fnMesh(dagPath);
//
//						// get the vertices that are part of the current mesh
//						MPointArray vertexList;
//						fnMesh.getPoints(vertexList, MSpace::kWorld);
//
//						// iterate through all the vertices
//						for (INT32 i = 0; i < vertexList.length(); i++)
//						{
//								vertexList[i].cartesianize();
//								MPoint point = vertexList[i];
//
//								// here is your data... now go do whatever you want with
//								// it. if you need a unique identifier for this vertex,
//								// use it's index in the mesh, and some kind of mesh id.
//								// these stay constant while exporting ( so long as the file is
//								// not edited )
//								(point.x, point.y, point.z);
//						}
//				}
//		}
//}
//
//// The creator method creates an instance of the circleNode class
//// and is the first method called by Maya
//// when a circleNode needs to be created.
////
//void* circle::creator()
//{
//		return new circle();
//}
//
//// The initialize routine is called after the node has been created.
//// It sets up the input and output attributes and adds them to the node.
//// Finally the dependencies are arranged so that when the inputs
//// change Maya knowns to call compute to recalculate the output values.
//// The inputs are: input, scale, frames
//// The outputs are: sineOutput, cosineOutput
////
//MStatus circle::initialize()
//{
//		MFnNumericAttribute nAttr;
//		MStatus				stat;
//
//		// Setup the input attributes
//		//
//		input = nAttr.create("input", "in", MFnNumericData::kFloat, 0.0,
//				&stat);
//		CHECK_MSTATUS(stat);
//		CHECK_MSTATUS(nAttr.setStorable(true));
//
//		scale = nAttr.create("scale", "sc", MFnNumericData::kFloat, 10.0,
//				&stat);
//		CHECK_MSTATUS(stat);
//		CHECK_MSTATUS(nAttr.setStorable(true));
//
//		frames = nAttr.create("frames", "fr", MFnNumericData::kFloat, 48.0,
//				&stat);
//		CHECK_MSTATUS(stat);
//		CHECK_MSTATUS(nAttr.setStorable(true));
//
//		// Setup the output attributes
//		//
//		sOutput = nAttr.create("sineOutput", "so", MFnNumericData::kFloat,
//				0.0, &stat);
//		CHECK_MSTATUS(stat);
//		CHECK_MSTATUS(nAttr.setWritable(false));
//		CHECK_MSTATUS(nAttr.setStorable(false));
//
//		cOutput = nAttr.create("cosineOutput", "co", MFnNumericData::kFloat,
//				0.0, &stat);
//		CHECK_MSTATUS(stat);
//		CHECK_MSTATUS(nAttr.setWritable(false));
//		CHECK_MSTATUS(nAttr.setStorable(false));
//
//		// Add the attributes to the node
//		//
//		CHECK_MSTATUS(addAttribute(input));
//		CHECK_MSTATUS(addAttribute(scale));
//		CHECK_MSTATUS(addAttribute(frames));
//		CHECK_MSTATUS(addAttribute(sOutput));
//		CHECK_MSTATUS(addAttribute(cOutput));
//
//		// Set the attribute dependencies
//		//
//		CHECK_MSTATUS(attributeAffects(input, sOutput));
//		CHECK_MSTATUS(attributeAffects(input, cOutput));
//		CHECK_MSTATUS(attributeAffects(scale, sOutput));
//		CHECK_MSTATUS(attributeAffects(scale, cOutput));
//		CHECK_MSTATUS(attributeAffects(frames, sOutput));
//		CHECK_MSTATUS(attributeAffects(frames, cOutput));
//		//extractVertices();
//		return MS::kSuccess;
//}
//
//// The constructor does nothing
////
//circle::circle() {}
//
//// The destructor does nothing
////
//circle::~circle() {}
//
//// The compute method is called by Maya when the input values
//// change and the output values need to be recomputed.
//// The input values are read then using sinf() and cosf()
//// the output values are stored on the output plugs.
////
//MStatus circle::compute(const MPlug& plug, MDataBlock& data)
//{
//
//		MStatus returnStatus;
//
//		// Check that the requested recompute is one of the output values
//		//
//		if (plug == sOutput || plug == cOutput) {
//				// Read the input values
//				//
//				MDataHandle inputData = data.inputValue(input, &returnStatus);
//				CHECK_MSTATUS(returnStatus);
//				MDataHandle scaleData = data.inputValue(scale, &returnStatus);
//				CHECK_MSTATUS(returnStatus);
//				MDataHandle framesData = data.inputValue(frames, &returnStatus);
//				CHECK_MSTATUS(returnStatus);
//
//				// Compute the output values
//				//
//				float currentFrame = inputData.asFloat();
//				float scaleFactor = scaleData.asFloat();
//				float framesPerCircle = framesData.asFloat();
//				float angle = 6.2831853f * (currentFrame / framesPerCircle);
//				float sinResult = sinf(angle) * scaleFactor;
//				float cosResult = cosf(angle) * scaleFactor;
//
//				// Store them on the output plugs
//				//
//				MDataHandle sinHandle = data.outputValue(circle::sOutput,
//						&returnStatus);
//				CHECK_MSTATUS(returnStatus);
//				MDataHandle cosHandle = data.outputValue(circle::cOutput,
//						&returnStatus);
//				CHECK_MSTATUS(returnStatus);
//				sinHandle.set(sinResult);
//				cosHandle.set(cosResult);
//				CHECK_MSTATUS(data.setClean(plug));
//		}
//		else {
//				return MS::kUnknownParameter;
//		}
//
//		return MS::kSuccess;
//}
//
//// The initializePlugin method is called by Maya when the circleNode
//// plugin is loaded.  It registers the circleNode which provides
//// Maya with the creator and initialize methods to be called when
//// a circleNode is created.
////
//MStatus initializePlugin(MObject obj)
//{
//		MStatus   status;
//		MFnPlugin plugin(obj, PLUGIN_COMPANY, "4.5", "Any");
//
//		status = plugin.registerNode("circle", circle::id,
//				circle::creator, circle::initialize);
//
//	//	status = plugin.register
//		if (!status) {
//				status.perror("registerNode");
//				return status;
//		}
//
//		return status;
//}
//
//// The unitializePlugin is called when Maya needs to unload the plugin.
//// It basically does the opposite of initialize by calling
//// the deregisterCommand to remove it.
////
//MStatus uninitializePlugin(MObject obj)
//{
//		MStatus   status;
//		MFnPlugin plugin(obj);
//
//		status = plugin.deregisterNode(circle::id);
//		if (!status) {
//				status.perror("deregisterNode");
//				return status;
//		}
//
//		return status;
//}
