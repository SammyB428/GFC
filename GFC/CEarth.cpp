/*
Author: Samuel R. Blackburn
Internet: wfc@pobox.com

"You can get credit for something or get it done, but not both."
Dr. Richard Garwin

The MIT License (MIT)

Copyright (c) 1997-2022 Sam Blackburn

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "GFC.h"
#pragma hdrstop

GeodesyFoundationClasses::CEarth::CEarth(Ellipsoid const ellipsoid_identifier )
{
   m_Initialize();
   SetEllipsoid( ellipsoid_identifier );
}

GeodesyFoundationClasses::CEarth::CEarth( GeodesyFoundationClasses::CEarth const& source )
{
   m_Initialize();
   Copy( source );
}

void GeodesyFoundationClasses::CEarth::AddLineOfSightDistanceAndDirectionToCoordinate( GeodesyFoundationClasses::CPolarCoordinate const& point_1, double const distance, double const direction, GeodesyFoundationClasses::CPolarCoordinate& point_2, double const height_above_surface_of_point_2 ) noexcept
{
   // The method used here is to convert the straight (line-of-sight) distance to
   // a surface distance and then find out the position using the surface distance.
   // This is a translation of the MMDIST routine found in the FORWRD3D program at
   // ftp://ftp.ngs.noaa.gov/pub/pcsoft/for_inv.3d/source/forwrd3d.for

   // Many thanks to Peter Dana (pdana@mail.utexas.edu) for educating me
   // on the finer points of Geodesy, one of which was how to compute
   // "second eccentricity squared"

   double const polar_eccentricity_squared  = ( ( m_EquatorialRadiusInMeters * m_EquatorialRadiusInMeters ) - ( m_PolarRadiusInMeters * m_PolarRadiusInMeters ) ) / ( m_PolarRadiusInMeters * m_PolarRadiusInMeters );
   double const point_1_latitude_in_radians = GeodesyFoundationClasses::CMath::ConvertDegreesToRadians( point_1.GetUpDownAngleInDegrees() );
   double const direction_in_radians        = GeodesyFoundationClasses::CMath::ConvertDegreesToRadians( direction );
   double const cosine_of_point_1_latitude  = GeodesyFoundationClasses::CMath::Cosine( point_1_latitude_in_radians );
   double const cosine_of_latitude_squared  = cosine_of_point_1_latitude * cosine_of_point_1_latitude;
   double const cosine_of_direction_squared = GeodesyFoundationClasses::CMath::Cosine( direction_in_radians ) * GeodesyFoundationClasses::CMath::Cosine( direction_in_radians );
   double const c = ( m_EquatorialRadiusInMeters * m_EquatorialRadiusInMeters ) / m_PolarRadiusInMeters;
   double const n = c / GeodesyFoundationClasses::CMath::SquareRoot( 1.0 + ( polar_eccentricity_squared * cosine_of_latitude_squared ) );
   double const r = n / ( 1.0 + ( polar_eccentricity_squared * cosine_of_latitude_squared * cosine_of_direction_squared ) );
   double const height_above_surface_of_point_1{ point_1.GetDistanceFromSurfaceInMeters() };
   double const difference_in_height = height_above_surface_of_point_2 - height_above_surface_of_point_1;
   double const term_1 = ( distance * distance ) - ( difference_in_height * difference_in_height );
   double const term_2 = 1.0 + ( height_above_surface_of_point_1 / r );
   double const term_3 = 1.0 + ( height_above_surface_of_point_2 / r );
   double const distance_adjusted_for_differences_in_height = GeodesyFoundationClasses::CMath::SquareRoot( term_1 / ( term_2 * term_3 ) );

   // printf( "distance_adjusted_for_differences_in_height is %.11lf\n", distance_adjusted_for_differences_in_height );

   double const two_r{ 2.0 * r };
   double const surface_distance{ two_r * GeodesyFoundationClasses::CMath::ArcSine(distance_adjusted_for_differences_in_height / two_r) };

   // printf( "surface_distance is %.11lf\n", surface_distance );

   AddSurfaceDistanceAndDirectionToCoordinate( point_1, surface_distance, direction, point_2 );
}

void GeodesyFoundationClasses::CEarth::AddSurfaceDistanceAndDirectionToCoordinate( GeodesyFoundationClasses::CEarthCoordinate const& point_1, double const distance, double const direction, CPolarCoordinate& point_2 ) noexcept
{
   GeodesyFoundationClasses::CPolarCoordinate polar_point_1;

   Convert( point_1, polar_point_1 );

   AddSurfaceDistanceAndDirectionToCoordinate( polar_point_1, distance, direction, point_2 );
}

void GeodesyFoundationClasses::CEarth::AddSurfaceDistanceAndDirectionToCoordinate( GeodesyFoundationClasses::CEarthCoordinate const& point_1, double const distance, double const direction, GeodesyFoundationClasses::CEarthCoordinate& point_2 ) noexcept
{
   GeodesyFoundationClasses::CPolarCoordinate polar_point_1;

   Convert( point_1, polar_point_1 );

   GeodesyFoundationClasses::CPolarCoordinate polar_point_2;

   AddSurfaceDistanceAndDirectionToCoordinate( polar_point_1, distance, direction, polar_point_2 );

   Convert( polar_point_2, point_2 );
}

void GeodesyFoundationClasses::CEarth::AddSurfaceDistanceAndDirectionToCoordinate( GeodesyFoundationClasses::CPolarCoordinate const& point_1, double const distance, double const direction, GeodesyFoundationClasses::CEarthCoordinate& point_2 ) noexcept
{
   GeodesyFoundationClasses::CPolarCoordinate polar_coordinate;

   AddSurfaceDistanceAndDirectionToCoordinate( point_1, distance, direction, polar_coordinate );

   Convert( polar_coordinate, point_2 );
}

void GeodesyFoundationClasses::CEarth::AddSurfaceDistanceAndDirectionToCoordinate( GeodesyFoundationClasses::CPolarCoordinate const& point_1, double const distance, double const direction, GeodesyFoundationClasses::CPolarCoordinate& point_2 ) noexcept
{
   // This is a translation of the Fortran routine DIRCT1 found in the
   // FORWRD3D program at:
   // ftp://ftp.ngs.noaa.gov/pub/pcsoft/for_inv.3d/source/forwrd3d.for

   double cosine_of_y = 0.0;
   double cz          = 0.0;
   double e           = 0.0;
   double sine_of_y   = 0.0;
   double term_1      = 0.0;
   double term_2      = 0.0;
   double term_3      = 0.0;
   double y           = 0.0;

   double const direction_in_radians = GeodesyFoundationClasses::CMath::ConvertDegreesToRadians( direction );
   double const eps = 0.000000000000005;
   double const r = 1.0 - m_Flattening;
   double const point_1_latitude_in_radians = GeodesyFoundationClasses::CMath::ConvertDegreesToRadians( point_1.GetUpDownAngleInDegrees() );
   double const point_1_longitude_in_radians = GeodesyFoundationClasses::CMath::ConvertDegreesToRadians( point_1.GetLeftRightAngleInDegrees() );

   double tangent_u = ( r * GeodesyFoundationClasses::CMath::Sine( point_1_latitude_in_radians ) ) / GeodesyFoundationClasses::CMath::Cosine( point_1_latitude_in_radians );

   double const sine_of_direction = GeodesyFoundationClasses::CMath::Sine( direction_in_radians );
   double const cosine_of_direction = GeodesyFoundationClasses::CMath::Cosine( direction_in_radians );
   double heading_from_point_2_to_point_1_in_radians{ 0.0 };

   if ( cosine_of_direction != 0.0 )
   {
      heading_from_point_2_to_point_1_in_radians = GeodesyFoundationClasses::CMath::ArcTangentOfYOverX( tangent_u, cosine_of_direction ) * 2.0;
   }

   double const cu = 1.0 / GeodesyFoundationClasses::CMath::SquareRoot( ( tangent_u * tangent_u ) + 1.0 );
   double const su = tangent_u * cu;
   double const sa = cu * sine_of_direction;
   double const c2a = ( (-sa) * sa ) + 1.0;
   double x = GeodesyFoundationClasses::CMath::SquareRoot( ( ( ( 1.0 / r / r ) - 1.0 ) * c2a ) + 1.0 ) + 1.0;
   x = ( x - 2.0 ) / x;
   double c = 1.0 - x;
   c = ( ( ( x * x ) / 4.0 ) + 1.0 ) / c;
   double d = ( ( 0.375 * ( x * x ) ) -1.0 ) * x;

   tangent_u = distance / r / m_EquatorialRadiusInMeters / c;

   y = tangent_u;

   bool exit_loop{ false };

   while( exit_loop != true )
   {
      sine_of_y = GeodesyFoundationClasses::CMath::Sine( y );
      cosine_of_y = GeodesyFoundationClasses::CMath::Cosine( y );
      cz = GeodesyFoundationClasses::CMath::Cosine( heading_from_point_2_to_point_1_in_radians + y );
      e = ( cz * cz * 2.0 ) - 1.0;
      c = y;
      x = e * cosine_of_y;
      y = ( e + e ) - 1.0;

      term_1 = ( sine_of_y * sine_of_y * 4.0 ) - 3.0;
      term_2 = ( ( term_1 * y * cz * d ) / 6.0 ) + x;
      term_3 = ( ( term_2 * d ) / 4.0 ) - cz;
      y = ( term_3 * sine_of_y * d ) + tangent_u;

      if (GeodesyFoundationClasses::CMath::AbsoluteValue( y - c ) > eps )
      {
         exit_loop = false;
      }
      else
      {
         exit_loop = true;
      }
   }

   heading_from_point_2_to_point_1_in_radians = ( cu * cosine_of_y * cosine_of_direction ) - ( su * sine_of_y );

   c = r * GeodesyFoundationClasses::CMath::SquareRoot( ( sa * sa ) + ( heading_from_point_2_to_point_1_in_radians * heading_from_point_2_to_point_1_in_radians ) );
   d = ( su * cosine_of_y ) + ( cu * sine_of_y * cosine_of_direction );

   double const point_2_latitude_in_radians = GeodesyFoundationClasses::CMath::ArcTangentOfYOverX( d, c );

   c = ( cu * cosine_of_y ) - ( su * sine_of_y * cosine_of_direction );
   x = GeodesyFoundationClasses::CMath::ArcTangentOfYOverX( sine_of_y * sine_of_direction, c );
   c = ( ( ( ( ( -3.0 * c2a ) + 4.0 ) * m_Flattening ) + 4.0 ) * c2a * m_Flattening ) / 16.0;
   d = ( ( ( ( e * cosine_of_y * c ) + cz ) * sine_of_y * c ) + y ) * sa;

   double const point_2_longitude_in_radians = ( point_1_longitude_in_radians + x ) - ( ( 1.0 - c ) * d * m_Flattening );

   heading_from_point_2_to_point_1_in_radians = GeodesyFoundationClasses::CMath::ArcTangentOfYOverX( sa, heading_from_point_2_to_point_1_in_radians ) + CMath::Pi();

   point_2.SetUpDownAngleInDegrees(GeodesyFoundationClasses::CMath::ConvertRadiansToDegrees( point_2_latitude_in_radians ) );
   point_2.SetLeftRightAngleInDegrees(GeodesyFoundationClasses::CMath::ConvertRadiansToDegrees( point_2_longitude_in_radians ) );
}

void GeodesyFoundationClasses::CEarth::Convert( GeodesyFoundationClasses::CEarthCoordinate const& cartesian_coordinate, GeodesyFoundationClasses::CPolarCoordinate& polar_coordinate ) const noexcept
{
   // convert from cartesian to polar

   auto const equatorial_radius_times_eccentricity_squared{ m_EquatorialRadiusInMeters * m_EccentricitySquared };

   auto const p{ CMath::SquareRoot((cartesian_coordinate.GetXCoordinateInMeters() * cartesian_coordinate.GetXCoordinateInMeters()) +
                                   (cartesian_coordinate.GetYCoordinateInMeters() * cartesian_coordinate.GetYCoordinateInMeters())) };

   double z_coordinate{ cartesian_coordinate.GetZCoordinateInMeters() }; // for convienance
   double const one_minus_eccentricity_squared{ 1.0 - m_EccentricitySquared };

   double temp_latitude{ z_coordinate / p / one_minus_eccentricity_squared };

   double old_value{ 0.0 };
   double part_a{ 0.0 };
   double part_b{ 0.0 };
   double part_c{ 0.0 };

   uint32_t loop_index{ 0 };
   uint32_t maximum_number_of_tries{ 1024 };

   bool convergence_was_acheived{ false };

   while( convergence_was_acheived != true && loop_index < maximum_number_of_tries )
   {
      old_value = temp_latitude;

      part_a = one_minus_eccentricity_squared * temp_latitude * temp_latitude;
      part_b = equatorial_radius_times_eccentricity_squared / GeodesyFoundationClasses::CMath::SquareRoot( 1.0 + part_a );
      part_c = p - part_b;

      temp_latitude = z_coordinate / part_c;

      loop_index++;

      if (GeodesyFoundationClasses::CMath::AbsoluteValue( temp_latitude - old_value ) > 0.000000000000000000001 )
      {
         // Oh well, try again...
      }
      else
      {
         // YIPEE!! We've reached convergence!
         convergence_was_acheived = true;
      }
   }

   if ( convergence_was_acheived == true )
   {
      // Save the UpDown angle in degrees
      
       double const latitude_angle_in_radians{ GeodesyFoundationClasses::CMath::ArcTangent(temp_latitude) };

      polar_coordinate.SetUpDownAngleInDegrees(GeodesyFoundationClasses::CMath::ConvertRadiansToDegrees( latitude_angle_in_radians ) ); // Latitude

      double const sine_of_latitude_in_radians   = GeodesyFoundationClasses::CMath::Sine(   latitude_angle_in_radians );
      double const cosine_of_latitude_in_radians = GeodesyFoundationClasses::CMath::Cosine( latitude_angle_in_radians );
      double const longitude_in_radians = GeodesyFoundationClasses::CMath::ArcTangentOfYOverX( cartesian_coordinate.GetYCoordinateInMeters(), cartesian_coordinate.GetXCoordinateInMeters() );

      polar_coordinate.SetLeftRightAngleInDegrees(GeodesyFoundationClasses::CMath::ConvertRadiansToDegrees( longitude_in_radians ) ); // Longitude

      auto const w{ GeodesyFoundationClasses::CMath::SquareRoot(1.0 - (m_EccentricitySquared * sine_of_latitude_in_radians * sine_of_latitude_in_radians)) };

      double const distance_from_center_to_surface_of_the_ellipsoid{ m_EquatorialRadiusInMeters / w };

      double distance_from_surface{ 0.0 };

      if (GeodesyFoundationClasses::CMath::AbsoluteValue( latitude_angle_in_radians ) < 0.7854 )
      {
         distance_from_surface = ( p / cosine_of_latitude_in_radians ) - distance_from_center_to_surface_of_the_ellipsoid;
      }
      else
      {
         distance_from_surface = ( z_coordinate / sine_of_latitude_in_radians ) - distance_from_center_to_surface_of_the_ellipsoid + ( m_EccentricitySquared * distance_from_center_to_surface_of_the_ellipsoid );
      }

      polar_coordinate.SetDistanceFromSurfaceInMeters( distance_from_surface );
   }
   else
   {
      // Oh well, we gave it a shot..
      polar_coordinate.Set( 0.0, 0.0, 0.0 );
   }
}

void GeodesyFoundationClasses::CEarth::Convert( GeodesyFoundationClasses::CPolarCoordinate const& polar_coordinate, GeodesyFoundationClasses::CEarthCoordinate& cartesian_coordinate ) const noexcept
{
   // convert from polar to cartesian

   auto const up_down_radians{ GeodesyFoundationClasses::CMath::ConvertDegreesToRadians(polar_coordinate.GetUpDownAngleInDegrees()) };// latitude
   auto const left_right_radians{ GeodesyFoundationClasses::CMath::ConvertDegreesToRadians(polar_coordinate.GetLeftRightAngleInDegrees()) };// longitude angle
   auto const sine_of_up_down_radians{ GeodesyFoundationClasses::CMath::Sine(up_down_radians) };
   auto const cosine_of_left_right_radians{ GeodesyFoundationClasses::CMath::Cosine(left_right_radians) };// cosine_of_longitude
   auto const cosine_of_up_down_radians{ GeodesyFoundationClasses::CMath::Cosine(up_down_radians) }; // cosine_of_latitude

   // Now we need to calculate the distance from the center of the ellipsoid to the surface of the ellipsoid
   auto const w{ GeodesyFoundationClasses::CMath::SquareRoot(1.0 - (m_EccentricitySquared * sine_of_up_down_radians * sine_of_up_down_radians)) };

   auto const distance_from_center_to_surface_of_the_ellipsoid{ m_EquatorialRadiusInMeters / w };

   double coordinate{ (distance_from_center_to_surface_of_the_ellipsoid + polar_coordinate.GetDistanceFromSurfaceInMeters()) * cosine_of_up_down_radians * cosine_of_left_right_radians };
   cartesian_coordinate.SetXCoordinateInMeters( coordinate );

   coordinate = ( distance_from_center_to_surface_of_the_ellipsoid + polar_coordinate.GetDistanceFromSurfaceInMeters() ) * cosine_of_up_down_radians * CMath::Sine( left_right_radians );
   cartesian_coordinate.SetYCoordinateInMeters( coordinate );

   coordinate = ( distance_from_center_to_surface_of_the_ellipsoid * ( 1.0 - m_EccentricitySquared ) + polar_coordinate.GetDistanceFromSurfaceInMeters() ) * sine_of_up_down_radians;
   cartesian_coordinate.SetZCoordinateInMeters( coordinate );
}

void GeodesyFoundationClasses::CEarth::Copy( GeodesyFoundationClasses::CEarth const& source ) noexcept
{
   if ( &source != this )
   {
      m_PolarRadiusInMeters      = source.m_PolarRadiusInMeters;
      m_EquatorialRadiusInMeters = source.m_EquatorialRadiusInMeters;
      m_Flattening               = source.m_Flattening;
      m_EccentricitySquared      = source.m_EccentricitySquared;
      m_EllipsoidID              = source.m_EllipsoidID;
   }
}

double GeodesyFoundationClasses::CEarth::GetLineOfSightDistanceFromCourse( GeodesyFoundationClasses::CEarthCoordinate const& current_location, GeodesyFoundationClasses::CEarthCoordinate const& point_a, GeodesyFoundationClasses::CEarthCoordinate const& point_b ) const noexcept
{
   // This function tells you how far off course you are from a straight line between
   // point_a and point_b.

   /*
   ** References:
   ** I got the formula from:
   ** Engineering Mathematics Handbook
   ** Jan J. Tuma, Ph.D.
   ** McGraw-Hill Book Company
   ** 1970
   ** Library of Congress Catalog Number 78-101174
   ** page 19, (a) Oblique triangle
   **
   ** Teach Yourself Trigonometry
   ** P. Abbott, B.A.
   ** English Universities Press Ltd.
   ** 102 Newgate Street
   ** London, E.C.I
   ** Originally published 1940
   ** I used the 1964 printing.
   ** Page 22, Figure 12 calls this "the altitude from the vertex A"
   */ 

   auto const distance_from_current_location_to_point_a{ GetLineOfSightDistance(current_location, point_a) };
   auto const distance_from_current_location_to_point_b{ GetLineOfSightDistance(current_location, point_b) };
   auto const distance_from_point_a_to_point_b{ GetLineOfSightDistance(point_a, point_b) };

   double p{ distance_from_current_location_to_point_a };
   p += distance_from_current_location_to_point_b;
   p += distance_from_point_a_to_point_b;
   p /= 2.0;

   double temp_double{ p };
   temp_double *= (double) ( p - distance_from_current_location_to_point_a );
   temp_double *= (double) ( p - distance_from_current_location_to_point_b );
   temp_double *= (double) ( p - distance_from_point_a_to_point_b );

   auto const area{ GeodesyFoundationClasses::CMath::SquareRoot(temp_double) };

   // The altitude from the vertex A is two times the area of the triangle divided by the baseline

   auto const distance_from_course = ( 2.0 * area ) / distance_from_point_a_to_point_b;

   return( distance_from_course );
}

double GeodesyFoundationClasses::CEarth::GetLineOfSightDistance( GeodesyFoundationClasses::CEarthCoordinate const& point_1, GeodesyFoundationClasses::CEarthCoordinate const& point_2 ) const noexcept
{
   // This function implements the Pythagoras method of computing the distance
   // between two points.
   // This is a line-of-sight algorithm. It does not take into account the 
   // curvature of the Earth. It is not a distance on the surface algorithm.
   // If you had a laser and connected the two points, this algorithm tells
   // you how long the laser beam is. 

   double x_coordinate{ point_1.GetXCoordinateInMeters() - point_2.GetXCoordinateInMeters() };
   double y_coordinate{ point_1.GetYCoordinateInMeters() - point_2.GetYCoordinateInMeters() };
   double z_coordinate{ point_1.GetZCoordinateInMeters() - point_2.GetZCoordinateInMeters() };

   // Square the coordinates
   x_coordinate *= x_coordinate;
   y_coordinate *= y_coordinate;
   z_coordinate *= z_coordinate;

   auto const distance = CMath::SquareRoot( x_coordinate + y_coordinate + z_coordinate );

   return( distance );
}

double GeodesyFoundationClasses::CEarth::GetLineOfSightDistance( GeodesyFoundationClasses::CEarthCoordinate const& point_1, GeodesyFoundationClasses::CPolarCoordinate const& point_2 ) const noexcept
{
   GeodesyFoundationClasses::CEarthCoordinate earth_center_earth_fixed_point_2;

   Convert( point_2, earth_center_earth_fixed_point_2 );

   return( GetLineOfSightDistance( point_1, earth_center_earth_fixed_point_2 ) );
}

double GeodesyFoundationClasses::CEarth::GetLineOfSightDistance( GeodesyFoundationClasses::CPolarCoordinate const& point_1, GeodesyFoundationClasses::CEarthCoordinate const& point_2 ) const noexcept
{
   GeodesyFoundationClasses::CEarthCoordinate earth_center_earth_fixed_point_1;

   Convert( point_1, earth_center_earth_fixed_point_1 );

   return( GetLineOfSightDistance( earth_center_earth_fixed_point_1, point_2 ) );
}

double GeodesyFoundationClasses::CEarth::GetLineOfSightDistance( GeodesyFoundationClasses::CPolarCoordinate const& point_1, GeodesyFoundationClasses::CPolarCoordinate const& point_2 ) const noexcept
{
   GeodesyFoundationClasses::CEarthCoordinate earth_center_earth_fixed_point_1;
   GeodesyFoundationClasses::CEarthCoordinate earth_center_earth_fixed_point_2;

   Convert( point_1, earth_center_earth_fixed_point_1 );
   Convert( point_2, earth_center_earth_fixed_point_2 );

   return( GetLineOfSightDistance( earth_center_earth_fixed_point_1, earth_center_earth_fixed_point_2 ) );
}

double GeodesyFoundationClasses::CEarth::GetSurfaceDistance( GeodesyFoundationClasses::CEarthCoordinate const& point_1, GeodesyFoundationClasses::CEarthCoordinate const& point_2, double * heading_from_point_1_to_point_2_p, double * heading_from_point_2_to_point_1_p ) const noexcept
{
   GeodesyFoundationClasses::CPolarCoordinate polar_point_1;
   GeodesyFoundationClasses::CPolarCoordinate polar_point_2;

   Convert( point_1, polar_point_1 );
   Convert( point_2, polar_point_2 );

   return( GetSurfaceDistance( polar_point_1, polar_point_2, heading_from_point_1_to_point_2_p, heading_from_point_2_to_point_1_p ) );
}

double GeodesyFoundationClasses::CEarth::GetSurfaceDistance( GeodesyFoundationClasses::CEarthCoordinate const& point_1, GeodesyFoundationClasses::CPolarCoordinate const& point_2, double * heading_from_point_1_to_point_2_p, double * heading_from_point_2_to_point_1_p ) const noexcept
{
   GeodesyFoundationClasses::CPolarCoordinate polar_point_1;

   Convert( point_1, polar_point_1 );

   return( GetSurfaceDistance( polar_point_1, point_2, heading_from_point_1_to_point_2_p, heading_from_point_2_to_point_1_p ) );
}

double GeodesyFoundationClasses::CEarth::GetSurfaceDistance( GeodesyFoundationClasses::CPolarCoordinate const& point_1, GeodesyFoundationClasses::CEarthCoordinate const& point_2, double * heading_from_point_1_to_point_2_p, double * heading_from_point_2_to_point_1_p ) const noexcept
{
   GeodesyFoundationClasses::CPolarCoordinate polar_point_2;

   Convert( point_2, polar_point_2 );

   return( GetSurfaceDistance( point_1, polar_point_2, heading_from_point_1_to_point_2_p, heading_from_point_2_to_point_1_p ) );
}

double GeodesyFoundationClasses::CEarth::GetSurfaceDistance( GeodesyFoundationClasses::CPolarCoordinate const& point_1, GeodesyFoundationClasses::CPolarCoordinate const& point_2, double * heading_from_point_1_to_point_2_p, double * heading_from_point_2_to_point_1_p ) const noexcept
{
   // This is a translation of the Fortran routine INVER1 found in the
   // INVERS3D program at:
   // ftp://ftp.ngs.noaa.gov/pub/pcsoft/for_inv.3d/source/invers3d.for

   // Eugeny Toukh [jecat@mail.ru] found a bug in the code.
   // Given the input valus of point1 1 deg S, 180 deg E and
   // point2 0 deg N, 1 deg W, we go into an endless loop.
   // The values of d and x will begin to oscillate, so he suggested
   // testing for this condition (if the current value of x equals
   // the last value of d) and exit the loop accordingly.

   // The ton of variables used...

   double c           = 0.0;
   double c2a         = 0.0;
   double cosine_of_x = 0.0;
   double cy          = 0.0;
   double cz          = 0.0;
   double d           = 0.0;
   double e           = 0.0;
   double sa          = 0.0;
   double sine_of_x   = 0.0;
   double sy          = 0.0;
   double y           = 0.0;
   double last_value_of_d = 0.00;

   // UpDown    == Latitude
   // LeftRight == Longitude

   auto const point_1_latitude_in_radians{ GeodesyFoundationClasses::CMath::ConvertDegreesToRadians(point_1.GetUpDownAngleInDegrees()) };
   auto const point_1_longitude_in_radians{ GeodesyFoundationClasses::CMath::ConvertDegreesToRadians(point_1.GetLeftRightAngleInDegrees()) };
   auto const point_2_latitude_in_radians{ GeodesyFoundationClasses::CMath::ConvertDegreesToRadians(point_2.GetUpDownAngleInDegrees()) };
   auto const point_2_longitude_in_radians{ GeodesyFoundationClasses::CMath::ConvertDegreesToRadians(point_2.GetLeftRightAngleInDegrees()) };

   double const r_value = 1.0 - m_Flattening;
   double tangent_1 = ( r_value * GeodesyFoundationClasses::CMath::Sine( point_1_latitude_in_radians ) ) / GeodesyFoundationClasses::CMath::Cosine( point_1_latitude_in_radians );
   double tangent_2 = ( r_value * GeodesyFoundationClasses::CMath::Sine( point_2_latitude_in_radians ) ) / GeodesyFoundationClasses::CMath::Cosine( point_2_latitude_in_radians );
   double c_value_1 = 1.0 / GeodesyFoundationClasses::CMath::SquareRoot( ( tangent_1 * tangent_1 ) + 1.0 );
   double s_value_1 = c_value_1 * tangent_1;
   double c_value_2 = 1.0 / GeodesyFoundationClasses::CMath::SquareRoot( ( tangent_2 * tangent_2 ) + 1.0 );
   double s = c_value_1 * c_value_2;

   double heading_from_point_2_to_point_1 = s * tangent_2; // backward_azimuth
   double heading_from_point_1_to_point_2 = heading_from_point_2_to_point_1 * tangent_1;

   double x = point_2_longitude_in_radians - point_1_longitude_in_radians;

   bool exit_loop{ false };

   while( exit_loop != true )
   {
      sine_of_x   = GeodesyFoundationClasses::CMath::Sine( x );
      cosine_of_x = GeodesyFoundationClasses::CMath::Cosine( x );
      tangent_1 = c_value_2 * sine_of_x;
      tangent_2 = heading_from_point_2_to_point_1 - ( s_value_1 * c_value_2 * cosine_of_x );
      sy = GeodesyFoundationClasses::CMath::SquareRoot( ( tangent_1 * tangent_1 ) + ( tangent_2 * tangent_2 ) );
      cy = ( s * cosine_of_x ) + heading_from_point_1_to_point_2;
      y = GeodesyFoundationClasses::CMath::ArcTangentOfYOverX( sy, cy );

      // Thanks to John Werner (werner@tij.wb.xerox.com) for
      // finding a bug where sy could be zero. Here's his fix:

      if ( ( s * sine_of_x ) == 0.0 && ( sy == 0.0 ) )
      {
         sa = 1.0;
      }
      else
      {
         sa = ( s * sine_of_x ) / sy;
      }

      c2a = ( (-sa) * sa ) + 1.0;
      cz = heading_from_point_1_to_point_2 + heading_from_point_1_to_point_2;

      if ( c2a > 0.0 )
      {
         cz = ( (-cz) / c2a ) + cy;
      }

      e = ( cz * cz * 2.0 ) - 1.0;
      c = ( ( ( ( ( -3.0 * c2a ) + 4.0 ) * m_Flattening ) + 4.0 ) * c2a * m_Flattening ) / 16.0;
      d = x;
      x = ( ( ( ( e * cy * c ) + cz ) * sy * c ) + y ) * sa;
      x = ( ( 1.0 - c ) * x * m_Flattening ) + point_2_longitude_in_radians - point_1_longitude_in_radians;

      if (GeodesyFoundationClasses::CMath::AbsoluteValue( d - x ) > 0.00000000000000000000005 && x != last_value_of_d )
      {
         exit_loop = false;
         last_value_of_d = d;
      }
      else
      {
         exit_loop = true;
      }
   }

   heading_from_point_1_to_point_2 = GeodesyFoundationClasses::CMath::ArcTangentOfYOverX( tangent_1, tangent_2 );

   auto temp_decimal_degrees = GeodesyFoundationClasses::CMath::ConvertRadiansToDegrees( heading_from_point_1_to_point_2 );

   if ( temp_decimal_degrees < 0.0 )
   {
      temp_decimal_degrees += 360.0;
   }

   if ( heading_from_point_1_to_point_2_p != nullptr )
   {
      // The user passed us a pointer, don't trust it.
      // If you are using Visual C++ on Windows NT, the following
      // try/catch block will ensure you won't blow up when random
      // pointers are passed to you. If you are on a legacy operating
      // system like Unix, you are screwed.

      try
      {
         *heading_from_point_1_to_point_2_p = temp_decimal_degrees;
      }
      catch( ... )
      {
         // Do Nothing
      }
   }

   heading_from_point_2_to_point_1 = GeodesyFoundationClasses::CMath::ArcTangentOfYOverX( c_value_1 * sine_of_x, ( (heading_from_point_2_to_point_1 * cosine_of_x ) - ( s_value_1 * c_value_2 ) ) ) + CMath::Pi();

   temp_decimal_degrees = GeodesyFoundationClasses::CMath::ConvertRadiansToDegrees( heading_from_point_2_to_point_1 );

   if ( temp_decimal_degrees < 0 )
   {
      temp_decimal_degrees += 360.0;
   }

   if ( heading_from_point_2_to_point_1_p != nullptr)
   {
      // The user passed us a pointer, don't trust it.
      // If you are using Visual C++ on Windows NT, the following
      // try/catch block will ensure you won't blow up when random
      // pointers are passed to you. If you are on a legacy operating
      // system like Unix, you are screwed.

      try
      {
         *heading_from_point_2_to_point_1_p = temp_decimal_degrees;
      }
      catch( ... )
      {
         // Do Nothing
      }
   }

   x = GeodesyFoundationClasses::CMath::SquareRoot( ( ( ( 1.0 / r_value / r_value ) - 1 ) * c2a ) + 1.0 ) + 1.0;
   x = ( x - 2.0 ) / x;
   c = 1.0 - x;
   c = ( ( ( x * x ) / 4.0 ) + 1.0 ) / c;
   d = ( ( 0.375 * ( x * x ) ) - 1.0 ) * x;

   // 1998-09-01
   // Thanks go to Gerard Murphy (bjmillar@dera.gov.uk) for finding a typo here.

   x = e * cy;

   s = ( 1.0 - e ) - e;

   double const term_1 = ( sy * sy * 4.0 ) - 3.0;
   double const term_2 = ( ( s * cz * d ) / 6.0 ) - x;
   double const term_3 = term_1 * term_2;
   double const term_4 = ( ( term_3 * d ) / 4.0 ) + cz;
   double const term_5 = ( term_4 * sy * d ) + y;

   s = term_5 * c * m_EquatorialRadiusInMeters * r_value;

   return( s );
}

void GeodesyFoundationClasses::CEarth::m_ComputeEccentricitySquared( void ) noexcept
{
   if ( m_Flattening == 0.0 )
   {
      m_EccentricitySquared = 0.0;
      return;
   }

   m_EccentricitySquared = ( 2.0 * m_Flattening ) - ( m_Flattening * m_Flattening );
}

void GeodesyFoundationClasses::CEarth::m_ComputeFlattening( void ) noexcept
{
   if ( m_EquatorialRadiusInMeters == 0.0 || m_PolarRadiusInMeters == 0.0 )
   {
      return;
   }

   m_Flattening = GeodesyFoundationClasses::CMath::AbsoluteValue( m_EquatorialRadiusInMeters - m_PolarRadiusInMeters ) / m_EquatorialRadiusInMeters;
}

void GeodesyFoundationClasses::CEarth::m_Initialize( void ) noexcept
{
   m_EllipsoidID              = Ellipsoid::Unknown;
   m_PolarRadiusInMeters      = 0.0;
   m_EquatorialRadiusInMeters = 0.0;
   m_Flattening               = 0.0;
   m_EccentricitySquared      = 0.0;
}

void GeodesyFoundationClasses::CEarth::SetEllipsoid(Ellipsoid const ellipsoid_identifier) noexcept
{
   m_EllipsoidID = ellipsoid_identifier;

   switch( ellipsoid_identifier )
   {
       case Ellipsoid::Perfect_Sphere:

         m_EquatorialRadiusInMeters = 6378137.0;
         m_PolarRadiusInMeters      = 6378137.0;

         break;

       case Ellipsoid::Airy:

         m_EquatorialRadiusInMeters = 6377563.396;
         m_PolarRadiusInMeters      = 6356256.909237;

         break;

       case Ellipsoid::Austrailian_National:

         m_EquatorialRadiusInMeters = 6378160.0;
         m_PolarRadiusInMeters      = 6356774.719195;

         break;

       case Ellipsoid::Bessell_1841:

         m_EquatorialRadiusInMeters = 6377397.155;
         m_PolarRadiusInMeters      = 6356078.962818;

         break;

       case Ellipsoid::Bessel_1841_Nambia:

         m_EquatorialRadiusInMeters = 6377483.865;
         m_PolarRadiusInMeters      = 6356165.382966;

         break;

       case Ellipsoid::Clarke_1866:

         m_EquatorialRadiusInMeters = 6378206.4;
         m_PolarRadiusInMeters      = 6356583.799999;

         break;

       case Ellipsoid::Clarke_1880:

         m_EquatorialRadiusInMeters = 6378249.145;
         m_PolarRadiusInMeters      = 6356514.86955;

         break;

       case Ellipsoid::Everest:

         m_EquatorialRadiusInMeters = 6377276.345;
         m_PolarRadiusInMeters      = 6356075.41314;

         break;

       case Ellipsoid::Fischer_1960_Mercury:

         m_EquatorialRadiusInMeters = 6378166.0;
         m_PolarRadiusInMeters      = 6356784.283607;

         break;

       case Ellipsoid::Fischer_1968:

         m_EquatorialRadiusInMeters = 6378150.0;
         m_PolarRadiusInMeters      = 6356768.337244;

         break;

       case Ellipsoid::GRS_1967:

         m_EquatorialRadiusInMeters = 6378160.0;
         m_PolarRadiusInMeters      = 6356774.516091;

         break;

       case Ellipsoid::GRS_1980:

         m_EquatorialRadiusInMeters = 6378137.0;
         m_PolarRadiusInMeters      = 6356752.31414;

         break;

       case Ellipsoid::Helmert_1906:

         m_EquatorialRadiusInMeters = 6378200.0;
         m_PolarRadiusInMeters      = 6356818.169628;

         break;

       case Ellipsoid::Hough:

         m_EquatorialRadiusInMeters = 6378270.0;
         m_PolarRadiusInMeters      = 6356794.343434;

         break;

       case Ellipsoid::International:

         m_EquatorialRadiusInMeters = 6378388.0;
         m_PolarRadiusInMeters      = 6356911.946128;

         break;

       case Ellipsoid::Krassovsky:

         m_EquatorialRadiusInMeters = 6378245.0;
         m_PolarRadiusInMeters      = 6356863.018773;

         break;

       case Ellipsoid::Modified_Airy:

         m_EquatorialRadiusInMeters = 6377340.189;
         m_PolarRadiusInMeters      = 6356034.447939;

         break;

       case Ellipsoid::Modified_Everest:

         m_EquatorialRadiusInMeters = 6377304.063;
         m_PolarRadiusInMeters      = 6356103.038993;

         break;

       case Ellipsoid::Modified_Fischer_1960:

         m_EquatorialRadiusInMeters = 6378155.0;
         m_PolarRadiusInMeters      = 6356773.320483;

         break;

       case Ellipsoid::South_American_1969:

         m_EquatorialRadiusInMeters = 6378160.0;
         m_PolarRadiusInMeters      = 6356774.719195;

         break;

       case Ellipsoid::Topex_Poseidon_Pathfinder_ITRF: // Source is http://neptune.gsfc.nasa.gov/~krachlin/corr/refframe.html

         m_EquatorialRadiusInMeters = 6378136.3;
         m_PolarRadiusInMeters      = 6356751.6005629376;

         break;

       case Ellipsoid::WGS_60:

         m_EquatorialRadiusInMeters = 6378165.0;
         m_PolarRadiusInMeters      = 6356783.286959;

         break;

       case Ellipsoid::WGS_66:

         m_EquatorialRadiusInMeters = 6378145.0;
         m_PolarRadiusInMeters      = 6356759.769489;

         break;

       case Ellipsoid::WGS_72:

         m_EquatorialRadiusInMeters = 6378135.0;
         m_PolarRadiusInMeters      = 6356750.520016;

         break;

       case Ellipsoid::WGS_84:

         // Computed polar radius from the flattening value specified at
         // http://acro.harvard.edu/SSA/BGA/wg84figs.html
         // because it had the most digits after the decimal point.

         m_EquatorialRadiusInMeters = 6378137.0;
         m_PolarRadiusInMeters      = 6356752.3142451793;

         break;

       case Ellipsoid::Unknown:
      default:

         m_EllipsoidID = Ellipsoid::Unknown;

         m_Initialize();
         return;
   }

   m_ComputeFlattening();
   m_ComputeEccentricitySquared();
}

void GeodesyFoundationClasses::CEarth::SetEllipsoidByRadii( double const equatorial_radius, double const polar_radius ) noexcept
{
   m_EquatorialRadiusInMeters = equatorial_radius;
   m_PolarRadiusInMeters      = polar_radius;
   m_EllipsoidID              = Ellipsoid::Custom;

   m_ComputeFlattening();
   m_ComputeEccentricitySquared();
}

void GeodesyFoundationClasses::CEarth::SetEllipsoidByEquatorialRadiusAndFlattening( double const equatorial_radius, double const flattening ) noexcept
{
   m_EquatorialRadiusInMeters = equatorial_radius;
   m_Flattening               = flattening;
   m_EllipsoidID              = Ellipsoid::Custom;

   // We must compute the polar radius

   double const temp_double = m_Flattening * m_EquatorialRadiusInMeters;

   m_PolarRadiusInMeters = m_EquatorialRadiusInMeters - temp_double;

   m_ComputeEccentricitySquared();
}

GeodesyFoundationClasses::CEarth& GeodesyFoundationClasses::CEarth::operator=( GeodesyFoundationClasses::CEarth const& source )
{
   Copy( source );
   return( *this );
}

#if 0
<HTML>
<HEAD><TITLE>GFC - CEarth</TITLE></HEAD>
<BODY>
<H1>CEarth</H1>
$Revision: 14 $
<HR>
<H2>Description</H2>
This class encapsulates the Earth. It holds the data necessary to perform the
calculations of distance and direction. The Earth is not a perfect sphere.
It is an ellipsoid (flattened at the top and bottom). All angles are expressed
in degrees and all distances are expressed in meters.
<H2>Methods</H2>
<DL COMPACT>
<DT><PRE>void <B>AddLineOfSightDistanceAndDirectionToCoordinate</B>( const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_1, double distance, double direction, <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_2, double height_above_surface_of_point_2 = 0.0 )</PRE><DD>
If you were to shine a laser from <CODE>point_1</CODE> pointing towards <CODE>direction</CODE>
for <CODE>distance</CODE> meters, this function will tell you what location you would
be at. It will also let you specify a point <CODE>height_above_surface_of_point_2</CODE>
meters above the surface. This method does not take into account the curvature of
the Earth.
<DT><PRE>void <B>AddSurfaceDistanceAndDirectionToCoordinate</B>( const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_1, double distance, double direction, <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_2 )
void <B>AddSurfaceDistanceAndDirectionToCoordinate</B>( const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_1, double distance, double direction, <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_2 )
void <B>AddSurfaceDistanceAndDirectionToCoordinate</B>( const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_1, double distance, double direction, <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_2 )
void <B>AddSurfaceDistanceAndDirectionToCoordinate</B>( const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_1, double distance, double direction, <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_2 )</PRE><DD>
This allows you to add a distance over the surface of the Earth to a location and
get a new location. It answers the question &quot;If I head out in this direction
for that amout of meters, where will I be?&quot;
<DT><PRE>void <B>Convert</B>( const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; cartesian_coordinate, <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; polar_coordinate ) const
void <B>Convert</B>( const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; polar_coordinate, <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; cartesian_coordinate ) const</PRE><DD>
This method allows you to convert from polar and cartestian coordinates.
<DT><PRE>double <B>GetDistanceToHorizon</B>( const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_1 ) const
double <B>GetDistanceToHorizon</B>( const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_1 ) const</PRE><DD>
This tells you how far (in meters) from the horizon <CODE>point_1</CODE> is.
<DT><PRE>double <B>GetEquatorialRadiusInMeters</B>( void ) const</PRE><DD>
This tells you what the equatorial radius is in meters for the selected ellipsoid.
<DT><PRE>double <B>GetPolarRadiusInMeters</B>( void ) const</PRE><DD>
This tells you what the polar radius is in meters for the selected ellipsoid.
<DT><PRE>double <B>GetLineOfSightDistanceFromCourse</B>( const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; current_location, const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_a, const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_b ) const</PRE><DD>
Draw a line from <CODE>point_a</CODE> to <CODE>point_b</CODE>. This function will tell
you how far <CODE>current_location</CODE> is from that line.
<DT><PRE>double <B>GetLineOfSightDistance</B>( const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_1, const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_2 ) const
double <B>GetLineOfSightDistance</B>( const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_1, const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_2 ) const
double <B>GetLineOfSightDistance</B>( const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_1, const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_2 ) const
double <B>GetLineOfSightDistance</B>( const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_1, const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_2 ) const</PRE><DD>
This will tell you how many meters it is between two points. It answers the question, &quot;If I
pointed a laser from <CODE>point_1</CODE> to <CODE>point_2</CODE>, how far would the
laser beam travel?&quot;
<DT><PRE>double <B>GetSurfaceDistance</B>( const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_1, const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_2, double * heading_from_point_1_to_point_2 = 0, double * heading_from_point_2_to_point_1 = 0 ) const
double <B>GetSurfaceDistance</B>( const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_1, const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_2, double * heading_from_point_1_to_point_2 = 0, double * heading_from_point_2_to_point_1 = 0 ) const
double <B>GetSurfaceDistance</B>( const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_1, const <A HREF="CEarthCoordinate.htm">CEarthCoordinate</A>&amp; point_2, double * heading_from_point_1_to_point_2 = 0, double * heading_from_point_2_to_point_1 = 0 ) const
double <B>GetSurfaceDistance</B>( const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_1, const <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A>&amp; point_2, double * heading_from_point_1_to_point_2 = 0, double * heading_from_point_2_to_point_1 = 0 ) const</PRE><DD>
This will tell you how many meters it is between two points. It answers the question, &quot;If I
were to walk from <CODE>point_1</CODE> to <CODE>point_2</CODE>, how far would I walk?&quot;
<DT><PRE>void <B>SetEllipsoid</B>( int ellipsoid )</PRE><DD>
This allows you to set the ellipsoid used by <B>CEarth</B> in its calculations. The
default is <CODE>WGS84</CODE> which is generally accepted as being the closest
approximation of the Earth's ellipsoid. The <CODE>ellipsoid</CODE> parameter
may be one of the following:
<UL>
<LI>Perfect_Sphere
<LI>Airy
<LI>Austrailian_National
<LI>Bessell_1841
<LI>Bessel_1841_Nambia
<LI>Clarke_1866
<LI>Clarke_1880
<LI>Everest
<LI>Fischer_1960_Mercury
<LI>Fischer_1968
<LI>GRS_1967
<LI>GRS_1980
<LI>Helmert_1906
<LI>Hough
<LI>International
<LI>Krassovsky
<LI>Modified_Airy
<LI>Modified_Everest
<LI>Modified_Fischer_1960
<LI>South_American_1969
<LI>Topex_Poseidon_Pathfinder_ITRF
<LI>WGS_60
<LI>WGS_66
<LI>WGS_72
<LI>WGS_84
<LI>NAD_27
<LI>Tokyo
</UL>
<DT><PRE>void <B>SetEllipsoidByRadii</B>( double equatorial_radius, double polar_radius )</PRE><DD>
This let's you use your own (custom) values to describe the ellipsoid of the Earth.
<DT><PRE>void <B>SetEllipsoidByEquatorialRadiusAndFlattening</B>( double equatorial_radius, double flattening )</PRE><DD>
This let's you use your own (custom) values to describe the ellipsoid of the Earth.
</DL>
<H2>Example</H2>
<PRE><CODE>#include &lt;stdio.h&gt;
#include &lt;GFC.h&gt;
#pragma hdrstop

void main( void )
{
   // Let's figure out how far it is from here to there

   <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A> here;
   <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A> there;

   // Convert from Latitude/Longitude to coordinates our system understands

   // here is 39 degrees 12.152 minutes North Latitude, 76 degrees 46.795 minutes West Longitude
   here.SetUpDownAngleInDegrees(     CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees(  39.0, 12.152, 0.0 ) );
   here.SetLeftRightAngleInDegrees(  CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees( -76.0, 46.795, 0.0 ) );

   // there is 12 degrees 8.535 minutes North Latitude, 68 degrees 16.547 West Longitude
   there.SetUpDownAngleInDegrees(    CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees(  12.0,  8.535, 0.0 ) );
   there.SetLeftRightAngleInDegrees( CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees( -68.0, 16.547, 0.0 ) );

   <B>CEarth</B> earth; // We are talking about the earth...

   double distance_in_meters         = 0.0;
   double heading_from_here_to_there = 0.0;
   double heading_from_there_to_here = 0.0;

   distance_in_meters = earth.GetSurfaceDistance( here, there, &heading_from_here_to_there, &heading_from_there_to_here );

   printf( &quot;Distance between here and there: %.23lf meters\nHeading from here to there:      %.19lf degrees\nHeading from there to here:      %.19lf degrees\n&quot;,
           distance_in_meters,
           heading_from_here_to_there,
           heading_from_there_to_here );

   double degrees = 0.0;
   double minutes = 0.0;
   double seconds = 0.0;

   CMath::ConvertDecimalDegreesToDegreesMinutesSeconds( heading_from_here_to_there, degrees, minutes, seconds );
   printf( &quot;Heading %lf degrees, %lf minutes, %lf seconds\n&quot;, degrees, minutes, seconds );

   CMath::ConvertDecimalDegreesToDegreesMinutesSeconds( heading_from_there_to_here, degrees, minutes, seconds );
   printf( &quot;Heading %lf degrees, %lf minutes, %lf seconds\n&quot;, degrees, minutes, seconds );
}</CODE></PRE>
<HR><I>Copyright, 1998, <A HREF="mailto:wfc@pobox.com">Samuel R. Blackburn</A></I><BR>
$Workfile: CEarth.cpp $<BR>
$Modtime: 8/06/00 10:30a $
</BODY>
</HTML>
ToolTipFormatLine=CEarth=m_EllipsoidID=<m_EllipsoidID>
#endif
