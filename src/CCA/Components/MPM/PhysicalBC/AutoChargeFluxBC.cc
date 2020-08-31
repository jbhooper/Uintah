/*
 * AutoChargeFluxBC.cc
 *
 *  Created on: Aug 28, 2020
 *      Author: jbhooper
 */

#include <CCA/Components/MPM/PhysicalBC/AutoChargeFluxBC.h>

#include <Core/Exceptions/ParameterNotFound.h>

#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>



namespace Uintah {

AutoChargeFluxBC::AutoChargeFluxBC(			ProblemSpecP	&	ps
								  ,const	GridP			&	grid
								  ,const	MPMFlags		*	f_MPMFlags
								  )
	{
		// Initialize everything to start
		d_surface = nullptr;
		d_loadCurve = nullptr;

		ProblemSpecP adult = ps->findBlock("geom_object");
		ProblemSpecP child = adult->findBlock();

		std::string go_type = child->getNodeName();
		d_surfaceType = go_type;

		switch (d_surfaceType) {
		case "box":
			d_surface = scinew BoxGeometryPiece(child);
			break;
		case "sphere":
			d_surface = scinew SphereGeometryPiece(child);
			break;
		case "cylinder":
			d_surface = scinew CylinderGeometryPiece(child);
			CylinderGeometryPiece* cgp = dynamic_cast<CylinderGeometryPiece*> (d_surface);
			if (cgp->axisymmetric_end()) {
				ps->require("res",d_res);
				Vector dx = grid->getLevel(0)->dCell();
				d_matlPointSize = Vector(dx.x()/((double) d_res.x()),
										 dx.y()/((double) d_res.y()),
										 dx.z()/((double) d_res.z()));
			}
			break;
		default:
			throw ParameterNotFound("* ERROR *: Surface geometry '" + go_type + "' not defined for AutoChargeFluxBC.",
									__FILE__, __LINE__);
		} // Geometry case switch

		d_numberMaterialPoints = 0;
		ps->get("numberOfParticlesOnLoadSurface", d_numberMaterialPoints);

		// Read and save th eload curve information
		d_loadCurve = scinew LoadCurve<double>(ps);

		// Warn if user's geometry object exceeds domain (there are reasons they might want to do this)
		Box boundingBox = d_surface->getBoundingBox();
		BBox computationalDomain;
		grid->getSpatialRange(computationalDomain);

		Point boundMinimum = boundingBox.lower();
		Point boundMaximum = boundingBox.upper();
		Point domainMinimum = computationalDomain.min();
		Point domainMaximum = computationalDomain.max();
		if ( boundMinimum.x() < domainMinimum.x() ) || ( boundMinimum.y() < domainMinimum.y() ) || (boundMinimum.z() < domainMinimum.z() )

	}

} // Namespace Uintah




