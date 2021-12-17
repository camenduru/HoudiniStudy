
#include "Mesh_Skeletonization.h"

// This is an automatically generated header file based on theDsFile, below,
// to provide SOP_StarParms, an easy way to access parameter values from
// SOP_StarVerb::cook with the correct type.
#include "Mesh_Skeletonization.proto.h"

#include <GU/GU_Detail.h>
#include <GEO/GEO_PrimPoly.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Director.h>

#include <PRM/PRM_Include.h>
#include <PRM/PRM_TemplateBuilder.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_StringHolder.h>
#include <SYS/SYS_Math.h>
#include <limits.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <fstream>

//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/extract_mean_curvature_flow_skeleton.h>
//#include <CGAL/boost/graph/split_graph_into_polylines.h>

using namespace HDK_Mesh_Skeletonization;

/// This is the internal name of the SOP type.
/// It isn't allowed to be the same as any other SOP's type name.
const UT_StringHolder Mesh_Skeletonization::theSOPTypeName("CGAL_Mesh_Skeletonization"_sh);

/// newSopOperator is the hook that Houdini grabs from this dll
/// and invokes to register the SOP.  In this case, we add ourselves
/// to the specified operator table.
void newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        Mesh_Skeletonization::theSOPTypeName,   // Internal name
        "CGAL_Mesh_Skeletonization",            // UI name
        Mesh_Skeletonization::myConstructor,    // How to build the SOP
        Mesh_Skeletonization::buildTemplates(), // My parameters
        0,										// Min # of sources
        1,										// Max # of sources
        nullptr,								// Custom local variables (none)
        OP_FLAG_GENERATOR));					// Flag it as generator
}

/// This is a multi-line raw string specifying the parameter interface
/// for this SOP.
static const char *theDsFile = R"THEDSFILE(
{
    name        parameters
	parm {
			name    "shourcept"      // Internal parameter name
			label   "Shourcept" // Descriptive parameter name for user interface
			type    integer
			default { "0" }     // Default for this parameter on new nodes
			range   { 0! 10 }   // The value is prevented from going below 2 at all.
								// The UI slider goes up to 50, but the value can go higher.
			export  all         // This makes the parameter show up in the toolbox
								// above the viewport when it's in the node's state.
		}
}
)THEDSFILE";

PRM_Template* Mesh_Skeletonization::buildTemplates()
{
	static PRM_TemplateBuilder templ("Mesh_Skeletonization.C"_sh, theDsFile);
	return templ.templates();
}

class Mesh_SkeletonizationVerb : public SOP_NodeVerb
{
public:
	Mesh_SkeletonizationVerb() {}
    virtual ~Mesh_SkeletonizationVerb() {}

    virtual SOP_NodeParms *allocParms() const { return new Mesh_SkeletonizationParms(); }
    virtual UT_StringHolder name() const { return Mesh_Skeletonization::theSOPTypeName; }
    virtual CookMode cookMode(const SOP_NodeParms *parms) const { return COOK_GENERIC; }

    virtual void cook(const CookParms &cookparms) const;

	/// This static data member automatically registers
	/// this verb class at library load time.
    static const SOP_NodeVerb::Register<Mesh_SkeletonizationVerb> theVerb;
};

// The static member variable definition has to be outside the class definition.
// The declaration is inside the class.
const SOP_NodeVerb::Register<Mesh_SkeletonizationVerb> Mesh_SkeletonizationVerb::theVerb;

const SOP_NodeVerb * Mesh_Skeletonization::cookVerb() const 
{ 
    return Mesh_SkeletonizationVerb::theVerb.get();
}

OP_ERROR Mesh_Skeletonization::Mesh_Skeletonization::cookMySop(OP_Context& context)
{
	// We must lock our inputs before we try to access their geometry.
	// OP_AutoLockInputs will automatically unlock our inputs when we return.
	// NOTE: Don't call unlockInputs yourself when using this!
	OP_AutoLockInputs inputs(this);
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();

	// Duplicate input geometry
	duplicateSource(0, context);

	// Flag the SOP as being time dependent (i.e. cook on time changes)
	flags().setTimeDep(true);

	return cookMyselfAsVerb(context);
}

typedef CGAL::Simple_cartesian<double>                        Kernel;
typedef Kernel::Point_3                                       Point;
typedef CGAL::Surface_mesh<Point>                             Triangle_mesh;
typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;
typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;

/// This is the function that does the actual work.
void Mesh_SkeletonizationVerb::cook(const SOP_NodeVerb::CookParms &cookparms) const
{
	Triangle_mesh tm;

	// This gives easy access to all of the current parameter values
	auto &&sopparms = cookparms.parms<Mesh_SkeletonizationParms>();
	GU_Detail *detail = cookparms.gdh().gdpNC();

	int32 sourcePtnum = sopparms.getShourcept();
	if (sourcePtnum < 0)
	{
		sourcePtnum = 0;
	}

	//detail->addFloatTuple(GA_ATTRIB_POINT, "distance", 1);
	//GA_RWHandleF DistanceHandl(detail, GA_ATTRIB_POINT, "distance");

	// NOTE: If you are only interested in the P attribute, use gdp->getP(),
	//       or don't bother making a new GA_RWHandleV3, and just use
	//       gdp->getPos3(ptoff) and gdp->setPos3(ptoff, Pvalue).
	//       This just gives an example supplying an attribute name.
	std::vector<vertex_descriptor> indexmap;
	GA_RWHandleV3 Phandle(detail->findAttribute(GA_ATTRIB_POINT, "P"));
	GA_Offset ptoff;
	GA_FOR_ALL_PTOFF(detail, ptoff)
	{
		UT_Vector3 Pvalue = Phandle.get(ptoff);
		Point ptTemp = Point(Pvalue.x(), Pvalue.y(), Pvalue.z());
		vertex_descriptor u = tm.add_vertex(ptTemp);
		indexmap.push_back(u);
	}
	
	std::vector<int32> ptoffsetmap;
	GEO_Primitive* prim;
	GA_FOR_ALL_PRIMITIVES(detail, prim)
	{
		GA_Range ptRange = prim->getPointRange(false);
		if (ptRange.getEntries() != 3)
		{
			cookparms.sopAddWarning(SOP_MESSAGE, "Input mesh must be Triangle_mesh!");
			return;
		}
		int32 size = 0;
		for (GA_Iterator it(ptRange); !it.atEnd(); ++it)
		{
			GA_Offset    ptoffset = (*it);
			ptoffsetmap.push_back((int32)ptoffset);
		}
	}
	for (size_t i = 0; i < ptoffsetmap.size() / 3; i++)
	{
		tm.add_face(indexmap[ptoffsetmap[i * 3]], indexmap[ptoffsetmap[i * 3 + 1]], indexmap[ptoffsetmap[i * 3 + 2]]);
	}

	Skeleton skeleton;
	Skeletonization mcs(tm);
	// 1. Contract the mesh by mean curvature flow.
	mcs.contract_geometry();
	// 2. Collapse short edges and split bad triangles.
	mcs.collapse_edges();
	mcs.split_faces();
	// 3. Fix degenerate vertices.
	mcs.detect_degeneracies();
	// Perform the above three steps in one iteration.
	if (sourcePtnum > 0)
	{
		for (int32 i = 0; i < sourcePtnum; i++)
		{
			mcs.contract();
		}
	}
	else
	{
		mcs.contract();
		//Iteratively apply step 1 to 3 until convergence.
		mcs.contract_until_convergence();
	}
	//	Convert the contracted mesh into a curve skeleton and
	// get the correspondent surface points
	mcs.convert_to_skeleton(skeleton);
	//std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
	//std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
	detail->clearAndDestroy();

	for (Skeleton_edge e : CGAL::make_range(edges(skeleton)))
	{
		const Point& s = skeleton[source(e, skeleton)].point;
		const Point& t = skeleton[target(e, skeleton)].point;
		GEO_PrimPoly *poly = (GEO_PrimPoly *)detail->appendPrimitive(GA_PRIMPOLY);
		
		int ptNumber = detail->appendPointOffset();
		UT_Vector3 PvalueNew(s.hx(), s.hy(), s.hz());
		Phandle.set(ptNumber, PvalueNew);
		poly->appendVertex(ptNumber);
		detail->bumpDataIdsForAddOrRemove(true, true, true);
		
		int ptNumber1 = detail->appendPointOffset();
		UT_Vector3 PvalueNew1(t.hx(), t.hy(), t.hz());
		Phandle.set(ptNumber1, PvalueNew1);
		poly->appendVertex(ptNumber1);
		detail->bumpDataIdsForAddOrRemove(true, true, true);
	}

}
