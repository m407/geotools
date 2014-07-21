/*
 *    GeoTools - The Open Source Java GIS Toolkit
 *    http://geotools.org
 *
 *    (C) 2006-2008, Open Source Geospatial Foundation (OSGeo)
 *
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation;
 *    version 2.1 of the License.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    This package contains formulas from the PROJ package of USGS.
 *    USGS's work is fully acknowledged here. This derived work has
 *    been relicensed under LGPL with Frank Warmerdam's permission.
 */
package org.geotools.referencing.operation.projection;

import java.awt.geom.Point2D;
import java.util.Collection;
import javax.measure.unit.NonSI;
import org.opengis.parameter.GeneralParameterDescriptor;
import org.opengis.parameter.ParameterDescriptor;
import org.opengis.parameter.ParameterDescriptorGroup;
import org.opengis.parameter.ParameterNotFoundException;
import org.opengis.parameter.ParameterValueGroup;
import org.opengis.referencing.operation.MathTransform;
import org.geotools.metadata.iso.citation.Citations;
import org.geotools.referencing.NamedIdentifier;
import org.geotools.resources.i18n.ErrorKeys;

import static java.lang.Math.*;


/**
 * Azimuthal Equidistant (EPSG code none).
 * <p>
 * <b>References:</b>
 * <ul>
 *   <li> A. Annoni, C. Luzet, E.Gubler and J. Ihde - Map Projections for Europe</li>
 *   <li> John P. Snyder (Map Projections - A Working Manual,
 *        U.S. Geological Survey Professional Paper 1395)</li>
 * </ul>
 *
 * @see <A HREF="http://mathworld.wolfram.com/AzimuthalEquidistantProjection.html">Azimuthal Equidistant Projection on MathWorld</A>
 * @see <A HREF="http://www.remotesensing.org/geotiff/proj_list/azimuthal_equidistant.html">"Azimuthal Equidistant" on RemoteSensing.org</A>
 *
 * @since 2.5
 * @version $Id$
 *
 *
 * @source $URL$
 * @author Jerry Huxtable  (for original code in Proj4)
 * @author Andrey Kuvshinov
 */
public class AzimuthalEquidistant extends MapProjection {
    /** For cross-version compatibility. */
    private static final long serialVersitonUID = 1736914748730174962L;

    /** Maximum difference allowed when comparing real numbers. */
    private static final double EPSILON = 1E-7;

    /** Epsilon for the comparaison of small quantities. */
    private static final double FINE_EPSILON = 1E-10;

    /** Epsilon for the comparaison of latitudes. */
    private static final double EPSILON_LATITUDE = 1E-10;
    
    /** The projection mode. */
    static final int OBLIQUE=0, EQUATORIAL=1, NORTH_POLE=2, SOUTH_POLE=3;

    /** The projection mode for this particular instance. */
    final int mode;

    /** Constant parameters. */
    private final static double TOL = 1.e-8;
    
    private double N1;
    private double Mp;
    private double He;
    private double G;
    private final double sinphi0, cosphi0;

    
    /**
     * Constructs a new map projection from the supplied parameters.
     *
     * @param  parameters The parameter values in standard units.
     * @throws ParameterNotFoundException if a mandatory parameter is missing.
     */
    protected AzimuthalEquidistant(final ParameterValueGroup parameters)
            throws ParameterNotFoundException
    {
        // Fetch parameters
        super(parameters);
        final Collection<GeneralParameterDescriptor> expected = getParameterDescriptors().descriptors();
        latitudeOfOrigin = doubleValue(expected, Provider.LATITUDE_OF_CENTRE,  parameters);
        centralMeridian  = doubleValue(expected, Provider.LONGITUDE_OF_CENTRE, parameters);
        ensureLatitudeInRange (Provider.LATITUDE_OF_CENTRE,  latitudeOfOrigin, true);
        ensureLongitudeInRange(Provider.LONGITUDE_OF_CENTRE, centralMeridian,  true);
        /*
         * Detects the mode (oblique, etc.).
         */
        final double t = abs(latitudeOfOrigin);
        if (abs(t - PI/2) < EPSILON_LATITUDE) {
            mode = latitudeOfOrigin < 0.0 ? SOUTH_POLE : NORTH_POLE;
        } else if (abs(t) < EPSILON_LATITUDE) {
            mode = EQUATORIAL;
        } else {
            mode = OBLIQUE;
        }

        switch (mode) {
            case NORTH_POLE: {
                sinphi0 = 1.0;
		cosphi0 = 0.0;
                break;
            }
            case SOUTH_POLE: {
                sinphi0 = -1.0;
		cosphi0 = 0.0;
                break;
            }
            case EQUATORIAL: {
                sinphi0 = 0.0;
                cosphi0 = 1.0;
                break;
            }
            case OBLIQUE: {
                sinphi0 = sin(latitudeOfOrigin);
		cosphi0 = cos(latitudeOfOrigin);
                break;
            }
            default: {
                throw new AssertionError(mode);
            }
        }
        if(!isSpherical){
            switch (mode) {
                case NORTH_POLE:
                        Mp = mlfn(PI/2, 1.0, 0.0);
                        break;
                case SOUTH_POLE:
                        Mp = mlfn(-PI/2, -1.0, 0.0);
                        break;
                case EQUATORIAL:  // Fall through
                case OBLIQUE:
                        N1 = 1.0 / sqrt(1.0 - excentricitySquared * sinphi0 * sinphi0);
                        G = sinphi0 * (excentricity / sqrt(1 - excentricitySquared));
                        He = cosphi0 * (excentricity / sqrt(1 - excentricitySquared));
                        break;
                default: {
                    throw new AssertionError(mode);
                }
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    public ParameterDescriptorGroup getParameterDescriptors() {
        return Provider.PARAMETERS;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ParameterValueGroup getParameterValues() {
        final ParameterValueGroup values = super.getParameterValues();
        final Collection<GeneralParameterDescriptor> expected = getParameterDescriptors().descriptors();
        set(expected, Provider.LATITUDE_OF_CENTRE,  values, latitudeOfOrigin);
        set(expected, Provider.LONGITUDE_OF_CENTRE, values, centralMeridian);
        return values;
    }

    /**
     * Transforms the specified (<var>&lambda;</var>,<var>&phi;</var>) coordinates
     * (units in radians) and stores the result in {@code ptDst} (linear distance
     * on a unit sphere).
     */
    protected Point2D transformNormalized(final double lambda, final double phi, Point2D ptDst)
            throws ProjectionException
    {
       
        double  rho, s, H, H2, c, Az, t, ct, st, cA, sA;
        final double x, y;

        double coslam = cos(lambda);
        double cosphi = cos(phi);
        double sinphi = sin(phi);
        
        switch (mode) {
            case NORTH_POLE:
                    coslam = - coslam;
            case SOUTH_POLE:
                    x = (rho = abs(Mp - mlfn(phi, sinphi, cosphi))) *
                            sin(lambda);
                    y = rho * coslam;
                    break;
            case EQUATORIAL:  //Fall through
            case OBLIQUE:
                    if (abs(lambda) < EPSILON_LATITUDE && abs(phi - latitudeOfOrigin) < EPSILON_LATITUDE) {
                            x = y = 0.0;
                            break;
                    }
                    t = atan2((1 - excentricitySquared) * sinphi + excentricitySquared * N1 * sinphi0 *
                            Math.sqrt(1. - excentricitySquared * sinphi * sinphi), cosphi);
                    ct = cos(t); st = Math.sin(t);
                    Az = atan2(Math.sin(lambda) * ct, cosphi0 * st - sinphi0 * coslam * ct);
                    cA = cos(Az); sA = sin(Az);
                    s = asin( abs(sA) < TOL ?
                            (cosphi0 * st - sinphi0 * coslam * ct) / cA :
                            Math.sin(lambda) * ct / sA );
                    H = He * cA;
                    H2 = H * H;
                    c = N1 * s * (1. + s * s * (- H2 * (1. - H2)/6. +
                            s * ( G * H * (1. - 2. * H2 * H2) / 8. +
                            s * ((H2 * (4. - 7. * H2) - 3. * G * G * (1. - 7. * H2)) /
                            120. - s * G * H / 48.))));
                    x = c * sA;
                    y = c * cA;
                    break;
            default: {
                throw new AssertionError(mode);
            }
        }

        if (ptDst != null) {
            ptDst.setLocation(x,y);
            return ptDst;
        }
        return new Point2D.Double(x,y);
    }

    /**
     * Transforms the specified (<var>x</var>,<var>y</var>) coordinate
     * and stores the result in {@code ptDst}.
     */
    @Override
    @SuppressWarnings("fallthrough")
    protected Point2D inverseTransformNormalized(double x, double y, Point2D ptDst)
            throws ProjectionException
    {
        double lambda, phi;

        double c, Az, cosAz, A, B, D, E, F, psi, t;

        int i;

        c = sqrt(x*x+y*y); //distance
        
        if ( c < EPSILON_LATITUDE) {
                phi = latitudeOfOrigin;
                lambda = 0.0;
        } else if (mode == OBLIQUE || mode == EQUATORIAL) {
                cosAz = cos(Az = atan2(x, y));
                t = cosphi0 * cosAz;
                B = excentricitySquared * t / (1 - excentricitySquared);
                A = - B * t;
                B *= 3.0 * (1.0 - A) * sinphi0;
                D = c / N1;
                E = D * (1.0 - D * D * (A * (1. + A) / 6. + B * (1.0 + 3.0*A) * D / 24.0));
                F = 1.0 - E * E * (A / 2.0 + B * E / 6.0);
                psi = asin(sinphi0 * cos(E) + t * sin(E));
                lambda = asin(sin(Az) * sin(E) / cos(psi));
                if ((t = abs(psi)) < EPSILON_LATITUDE)
                        phi = 0.0;
                else if (abs(t - PI/2) < 0.0)
                        phi = PI/2;
                else
                        phi = atan((1.0 - excentricitySquared * F * sinphi0 / sin(psi)) * tan(psi) / (1 - excentricitySquared));
        } else {
                phi = inv_mlfn(mode == NORTH_POLE ? Mp - c : Mp + c);
                lambda = atan2(x, mode == NORTH_POLE ? -y : y);
        }

        
        if (ptDst != null) {
            ptDst.setLocation(lambda, phi);
            return ptDst;
        }
        return new Point2D.Double(lambda, phi);
    }


    /**
     * Provides the transform equations for the spherical case.
     *
     * @version $Id$
     * @author Martin Desruisseaux
     */
    private static final class Spherical extends AzimuthalEquidistant {
        /**
         * For cross-version compatibility.
         */
        private static final long serialVersionUID = 231163136379134540L;

        /**
         * Constructs a new map projection from the suplied parameters.
         *
         * @param  parameters The parameter values in standard units.
         * @throws ParameterNotFoundException if a mandatory parameter is missing.
         */
        protected Spherical(final ParameterValueGroup parameters)
                throws ParameterNotFoundException
        {
            super(parameters);
            ensureSpherical();
        }

        /**
         * Transforms the specified (<var>&lambda;</var>,<var>&phi;</var>) coordinates
         * (units in radians) and stores the result in {@code ptDst} (linear distance
         * on a unit sphere).
         */
        @Override
        protected Point2D transformNormalized(final double lambda, final double phi, Point2D ptDst)
                throws ProjectionException
        {
            // Compute using ellipsoidal formulas, for comparaison later.
            assert (ptDst = super.transformNormalized(lambda, phi, ptDst)) != null;

            double sinphi = sin(phi);
            double cosphi = cos(phi);
            double coslam = cos(lambda);
            double x,y;
            
            switch (mode) {
                case EQUATORIAL: //Fall throughEPS10
                case OBLIQUE:
                    y = sin(latitudeOfOrigin) * sinphi + cos(latitudeOfOrigin) * cosphi * coslam;
                    if (abs(abs(y) - 1.0) < TOL)
                            if (y < 0.0)
                                    throw new ProjectionException(); 
                            else
                                    x = y = 0.;
                    else {
                            y = acos(y);
                            y /= sin(y);
                            x = y * cosphi * sin(lambda);
                            y *=  cos(latitudeOfOrigin) * sinphi - sin(latitudeOfOrigin) * cosphi * coslam;
                    }
                    break;
                case NORTH_POLE:
                    if (abs(phi + PI/2) < EPSILON_LATITUDE)
                            throw new ProjectionException();
                    x = (y = (PI/2 - phi)) * sin(lambda);
                    y *= (-coslam);
                    break;
                case SOUTH_POLE:
                    if (abs(phi - PI/2) < EPSILON_LATITUDE)
                            throw new ProjectionException();
                    x = (y = (PI/2 + phi)) * sin(lambda);
                    y *= coslam;
                    break;
                default: {
                    throw new AssertionError(mode);
                }

            }

            assert checkTransform(x, y, ptDst);
            if (ptDst != null) {
                ptDst.setLocation(x,y);
                return ptDst;
            }
            return new Point2D.Double(x,y);
        }

        /**
         * Transforms the specified (<var>x</var>,<var>y</var>) coordinate
         * and stores the result in {@code ptDst} using equations for a sphere.
         */
        @Override
        protected Point2D inverseTransformNormalized(double x, double y, Point2D ptDst)
                throws ProjectionException
        {
            // Compute using ellipsoidal formulas, for comparaison later.
            assert (ptDst = super.inverseTransformNormalized(x, y, ptDst)) != null;

            double lambda, phi;
            
            double cosc, c_rh, sinc;
            
            c_rh = sqrt(x*x+y*y); //distance
            
            if ( c_rh > PI) {
                    if (c_rh - EPSILON_LATITUDE > PI)
                            throw new ProjectionException(); 
                    c_rh = PI;
            } else if (c_rh < EPSILON_LATITUDE) {
                    phi = latitudeOfOrigin;
                    lambda = 0.0;
                    
                    assert checkInverseTransform(lambda, phi, ptDst);
                    if (ptDst != null) {
                        ptDst.setLocation(lambda, phi);
                        return ptDst;
                    }
                    return new Point2D.Double(lambda, phi);
            } 
            
            if (mode == OBLIQUE || mode == EQUATORIAL) {
                    sinc = sin(c_rh);
                    cosc = cos(c_rh);
                    if (mode == EQUATORIAL) {
                            phi = asin(y * sinc / c_rh);
                            x *= sinc;
                            y = cosc * c_rh;
                    } else {
                            phi = asin(cosc * sin(latitudeOfOrigin) + y * sinc * cos(latitudeOfOrigin) /
                                    c_rh);
                            y = (cosc - sin(latitudeOfOrigin) * sin(phi)) * c_rh;
                            x *= sinc * cos(latitudeOfOrigin);
                    }
                    lambda = y == 0.0 ? 0.0 : atan2(x, y);
            } else if (mode == NORTH_POLE) {
                    phi = PI/2 - c_rh;
                    lambda = atan2(x, -y);
            } else {
                    phi = c_rh - PI/2;
                    lambda = atan2(x, y);
            }

            assert checkInverseTransform(lambda, phi, ptDst);
            if (ptDst != null) {
                ptDst.setLocation(lambda, phi);
                return ptDst;
            }
            return new Point2D.Double(lambda, phi);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    ////////                                                                          ////////
    ////////                                 PROVIDERS                                ////////
    ////////                                                                          ////////
    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    /**
     * The {@linkplain org.geotools.referencing.operation.MathTransformProvider math transform
     * provider} for an {@linkplain Azimuthal Equidistant} projection
     * (EPSG code none).
     *
     * @since 2.5
     * @version $Id$
     * @author Beate Stollberg
     *
     * @see org.geotools.referencing.operation.DefaultMathTransformFactory
     */
    public static class Provider extends AbstractProvider {
        /**
         * For cross-version compatibility.
         */
        private static final long serialVersionUID = 3877793025552244132L;

        /**
         * The operation parameter descriptor for the {@link #latitudeOfOrigin}
         * parameter value. Valid values range is from -90 to 90°. Default value is 0.
         */
        public static final ParameterDescriptor LATITUDE_OF_CENTRE = createDescriptor(
                new NamedIdentifier[] {
                    new NamedIdentifier(Citations.OGC,      "latitude_of_center"),
                    new NamedIdentifier(Citations.EPSG,     "Latitude of natural origin"),
                    new NamedIdentifier(Citations.EPSG,     "Spherical latitude of origin"),
                    new NamedIdentifier(Citations.ESRI,     "Latitude_Of_Origin"),
                    new NamedIdentifier(Citations.GEOTIFF,  "ProjCenterLat")
                },
                0, -90, 90, NonSI.DEGREE_ANGLE);

        /**
         * The operation parameter descriptor for the {@link #centralMeridian}
         * parameter value. Valid values range is from -180 to 180°. Default value is 0.
         */
        public static final ParameterDescriptor LONGITUDE_OF_CENTRE = createDescriptor(
                new NamedIdentifier[] {
                    new NamedIdentifier(Citations.OGC,      "longitude_of_center"),
                    new NamedIdentifier(Citations.EPSG,     "Longitude of natural origin"),
                    new NamedIdentifier(Citations.EPSG,     "Spherical longitude of origin"),
                    new NamedIdentifier(Citations.ESRI,     "Central_Meridian"),
                    new NamedIdentifier(Citations.GEOTIFF,  "ProjCenterLong")
                },
                0, -180, 180, NonSI.DEGREE_ANGLE);

        /**
         * The parameters group.
         */
        static final ParameterDescriptorGroup PARAMETERS = createDescriptorGroup(new NamedIdentifier[] {
            new NamedIdentifier(Citations.OGC,     "Azimuthal_Equidistant"),
            new NamedIdentifier(Citations.EPSG,    "Azimuthal Equidistant"),
            new NamedIdentifier(Citations.EPSG,    "Azimuthal Equidistant (Spherical)"),
            new NamedIdentifier(Citations.GEOTIFF, "CT_AzimuthalEquidistant"),
            //new NamedIdentifier(Citations.EPSG,    "9820"),
        },  new ParameterDescriptor[] {
                SEMI_MAJOR,         SEMI_MINOR,
                LATITUDE_OF_CENTRE, LONGITUDE_OF_CENTRE,
                FALSE_EASTING,      FALSE_NORTHING
        });

        /**
         * Constructs a new provider.
         */
        public Provider() {
            super(PARAMETERS);
        }

        /**
         * Creates a transform from the specified group of parameter values.
         *
         * @param  parameters The group of parameter values.
         * @return The created math transform.
         * @throws ParameterNotFoundException if a required parameter was not found.
         */
        public MathTransform createMathTransform(final ParameterValueGroup parameters)
                throws ParameterNotFoundException
        {
            return isSpherical(parameters) ? new Spherical(parameters) :
                    new AzimuthalEquidistant(parameters);
        }
    }
}
