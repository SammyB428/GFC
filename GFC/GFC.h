/*
Author: Samuel R. Blackburn
Internet: wfc@pobox.com

"You can get credit for something or get it done, but not both."
Dr. Richard Garwin

The MIT License (MIT)

Copyright (c) 1997-2023 Sam Blackburn

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

/* SPDX-License-Identifier: MIT */

// Definitions
// Ellipsoid
//  A flattened sphere. Take a basketball (a sphere), let some air out of it then
//  stand on it. You now have an ellipsoid. The basketball is flat on the bottom
//  and the top but still round in the middle.
// Equator
//  0 degrees Latitude
// Flattening
//  This is a ratio between polar radius and equatorial radius. Other names for
//  flattening are elliptocyte or oblateness.
// Latitude
//  The lines that run around the Earth like belts. They measure up/down angles.
// Longituide
//  The lines that slice the Earth like an orange from top to bottom. They
//  measure left/right angles.
// Meridian
//   One of the imaginary circles on the earth's surface passing through the north
//   and south poles 
// Pi
//  The most famous constant. It is a ratio.
//  It is the Circumference of a circle divided by
//  the diameter of the circle. It is 3.14 (or roughly 22 divivded by 7).
// Prime Meridian 
//  0 degrees Longitude
// Radian
//  The unit of measurement of a circle. It is a ratio.
//  It's kind of hard to explain what
//  a radian is so let me give you an example. Take the basketball and cut
//  a length of string equal to the radius (i.e. half the width) of the
//  ball. Now, lay this string on top of the ball. See how it forms an arc?
//  Now, draw a line from each end of the string to the center of the ball.
//  The angle of these two lines at the center of the ball is roughly
//  57.2958 degrees. This is true no matter what the size of the ball is.
//  Why? Because as the size of the ball increases, so does the radius.
//  So, a radian can be considered as a relationship between the radius of
//  a sphere (ball) and the circumference of the sphere. Now, there's something
//  interesting about the number of degrees in a radian (52.2958). It is
//  equal to 180 divided by Pi. Yup! Go Figure! OK. Now we're getting somewhere.
//  To find the length of the arc when all you have is the number of radians
//  in the angle at the center of the sphere, you multiply the radius times
//  the number of radians. It's that simple (on a perfect sphere anyway).

// Geodetic Datum - description of the shape of the earth

// Reference Ellipsoid - A flattened sphere. Let some air out of
// a basketball, then stand on it. The ball is now an ellipsoid.
// Your feet are standing on the North Pole and the floor is the
// South Pole.

// Ellipsoids are described by their polar and equatorial radii.
// Polar radius is also known as semi-minor axis.
// Equatorial radius is also know as semi-major axis.
// All other items of interest about the ellipsoid are derived
// from these two data.
// Flattening is ( Equatorial radius - Polar radius ) / Equatorial radius
// There's another thing called First Eccentricity Squared, this is computed
// by ( 2 * Flattening ) - ( Flattening squared ).

// Coordinates - a means of expressing a point in space
//
// Cartesian coordinates are X, Y, and Z coordinates. These are always
// expressed in distances from 0, 0, 0 (i,e, the center of the earth).
//
// Polar coordinates are theta (often seen as a letter O with a line
// through its middle running left-to-right), phi (often seen as a
// letter O with a line through its middle running top-to-bottom), and
// r (a distance).
// These are two angles and a distance. The angles are measured from
// a common origin (the center of the earth). Theta is the plane that
// runs through the equator, phi runs through the poles. R is the 
// distance along the line running along the phi angle. For simplicity
// sake, we will refer to theta as Equatorial angle, phi as the
// polar angle, and r as Polar Distance.
//
// Converting coordinates
//
// You can convert from polar to cartesian cordinates using the following
// formula:
// X = Polar distance * cosine of Polar angle * cosine of Equatorial angle
// Y = Polar distance * cosine of Polar angle * sine of Equatorial angle
// Z = Polar distance * sine of Polar angle

// Applying this to the real world
//
// Cartesian coordinates ar commonly refered to an ECEF X,Y,Z coordinates.
// This means Earth Centered, Earth Fixed. The 0,0,0 coordinate is the
// center of the earth. The X axis runs towards the Prime Meridian.
// The Y axis runs towards the equator and the Z axis runs toward
// the poles. Positive Z numbers go towards the North pole while
// negative numbers go towards the South Pole. Positive 

// Computing Distance
//
// If you have two cartesian coordinates, you can compute a line
// of sight (as the bullet flies, aiming a laser, pointing in a straight line)
// by this formula (some guy named Pythagoras figured this out):
// SQRT( ( X1 - X2 ) ^ 2 + ( Y1 - Y2 ) ^ 2 + ( Z1 - Z2 ) ^ 2 )
//
// or in pseudo code:
//
// cartesian_coordinate first_location;
// cartesian_coordinate second_location;
//
// double temp_x;
// double temp_y;
// double temp_z;
//
// temp_x = first_location.X - second_location.X;
// temp_y = first_location.Y - second_location.Y;
// temp_z = first_location.Z - second_location.Z;
//
// temp_x = temp_x * temp_x; // square them
// temp_y = temp_y * temp_y;
// temp_z = temp_z * temp_z;
//
// double temp_double;
//
// temp_double = temp_x + temp_y + temp_z;
//
// double distance;
//
// distance = sqrt( temp_double );
//
// While this will give you distance, it will not give you direction.

#include <cstdint>
#include <cmath>

namespace GeodesyFoundationClasses
{
    class CMath
    {
        // This class encapsulates all of the math functions. It is here to allow you
        // to replace the C Runtime functions with your home-grown (and maybe better
        // implementation) routines

        public:

        [[nodiscard]] static inline double AbsoluteValue(double const value) noexcept;
        [[nodiscard]] static inline double ArcCosine(double const value) noexcept;
        [[nodiscard]] static inline double ArcSine(double const value) noexcept;
        [[nodiscard]] static inline double ArcTangent(double const value) noexcept;
        [[nodiscard]] static inline double ArcTangentOfYOverX(double const y, double const x) noexcept;
        [[nodiscard]] static inline double Ceiling(double const value) noexcept;
        [[nodiscard]] static inline constexpr double ConvertDegreesToRadians(double const degrees) noexcept;
        [[nodiscard]] static inline constexpr double ConvertRadiansToDegrees(double const radians) noexcept;
        static inline void ConvertDecimalDegreesToDegreesMinutesSeconds(double const decimal_degrees, double& degrees, double& minutes, double& seconds) noexcept;
        // West is negative, East is positive, North is positive, south is negative
        [[nodiscard]] static inline constexpr double ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees(double const degrees, double const minutes, double const seconds) noexcept;
        [[nodiscard]] static inline double Cosine(double const value) noexcept;
        [[nodiscard]] static inline double HyperbolicCosine(double const value) noexcept;
        [[nodiscard]] static inline constexpr double Pi(void) noexcept;
        [[nodiscard]] static inline double Sine(double const value) noexcept;
        [[nodiscard]] static inline double SquareRoot(double const value) noexcept;
    };

#include "CMath.inl" // Implementations of the inline functions

    class CEarthCoordinate
    {
        // This is a Cartesian coordinate (Earth Centered, Earth Fixed)

        protected:

        double m_X_CoordinateInMeters{ 0.0 }; // Positive points to intersection of the Prime Meridian and the equator
        double m_Y_CoordinateInMeters{ 0.0 }; // Positive points to the intersection of 90 degrees east of Prime Meridian and the equator
        double m_Z_CoordinateInMeters{ 0.0 }; // Positive points to the North Pole, Negative towards the South Pole

        public:

        CEarthCoordinate() noexcept
        {
            m_X_CoordinateInMeters = 0.0;
            m_Y_CoordinateInMeters = 0.0;
            m_Z_CoordinateInMeters = 0.0;
        }

        CEarthCoordinate(CEarthCoordinate const& source) noexcept
        {
            Copy(source);
        }

        ~CEarthCoordinate() = default;

        inline constexpr void Copy(CEarthCoordinate const& source) noexcept
        {
            m_X_CoordinateInMeters = source.m_X_CoordinateInMeters;
            m_Y_CoordinateInMeters = source.m_Y_CoordinateInMeters;
            m_Z_CoordinateInMeters = source.m_Z_CoordinateInMeters;
        }

        inline constexpr void Get(double& x_coordinate, double& y_coordinate, double& z_coordinate) const noexcept
        {
            x_coordinate = m_X_CoordinateInMeters;
            y_coordinate = m_Y_CoordinateInMeters;
            z_coordinate = m_Z_CoordinateInMeters;
        }

        [[nodiscard]] inline constexpr double GetXCoordinateInMeters(void) const noexcept
        {
            return(m_X_CoordinateInMeters);
        }

        [[nodiscard]] inline constexpr double GetYCoordinateInMeters(void) const noexcept
        {
            return(m_Y_CoordinateInMeters);
        }

        [[nodiscard]] inline constexpr double GetZCoordinateInMeters(void) const noexcept
        {
            return(m_Z_CoordinateInMeters);
        }

        inline constexpr void SetXCoordinateInMeters(double const x_coordinate) noexcept
        {
            m_X_CoordinateInMeters = x_coordinate;
        }

        inline constexpr void SetYCoordinateInMeters(double const y_coordinate) noexcept
        {
            m_Y_CoordinateInMeters = y_coordinate;
        }

        inline constexpr void SetZCoordinateInMeters(double const z_coordinate) noexcept
        {
            m_Z_CoordinateInMeters = z_coordinate;
        }

        inline constexpr void Set(double const x_coordinate, double const y_coordinate, double const z_coordinate) noexcept
        {
            m_X_CoordinateInMeters = x_coordinate;
            m_Y_CoordinateInMeters = y_coordinate;
            m_Z_CoordinateInMeters = z_coordinate;
        }

        CEarthCoordinate& operator = (CEarthCoordinate const& source)
        {
            Copy(source);
            return(*this);
        }
    };

    class CPolarCoordinate
    {
        protected:

        double m_UpDownAngleInDegrees{ 0.0 }; // Polar Angle, Phi
        double m_LeftRightAngleInDegrees{ 0.0 }; // Equatorial Angle, Theta
        double m_DistanceFromSurfaceInMeters{ 0.0 };

        public:

        CPolarCoordinate() noexcept
        {
            m_UpDownAngleInDegrees = 0.0;
            m_LeftRightAngleInDegrees = 0.0;
            m_DistanceFromSurfaceInMeters = 0.0;
        }

        CPolarCoordinate(CPolarCoordinate const& source) noexcept
        {
            Copy(source);
        }

        ~CPolarCoordinate() = default;

        inline constexpr void Copy(CPolarCoordinate const& source) noexcept
        {
            m_UpDownAngleInDegrees = source.m_UpDownAngleInDegrees;
            m_LeftRightAngleInDegrees = source.m_LeftRightAngleInDegrees;
            m_DistanceFromSurfaceInMeters = source.m_DistanceFromSurfaceInMeters;
        }

        inline constexpr void Get(double& up_down_angle, double& left_right_angle, double& distance_from_surface) const noexcept
        {
            up_down_angle = m_UpDownAngleInDegrees;
            left_right_angle = m_LeftRightAngleInDegrees;
            distance_from_surface = m_DistanceFromSurfaceInMeters;
        }

        [[nodiscard]] inline constexpr double GetUpDownAngleInDegrees(void) const
        {
            return(m_UpDownAngleInDegrees);
        }

        [[nodiscard]] inline constexpr double GetLatitude(void) const
        {
            return(GetUpDownAngleInDegrees());
        }

        [[nodiscard]] inline constexpr double GetLeftRightAngleInDegrees(void) const noexcept
        {
            return(m_LeftRightAngleInDegrees);
        }

        [[nodiscard]] inline constexpr double GetLongitude(void) const
        {
            return(GetLeftRightAngleInDegrees());
        }

        [[nodiscard]] inline constexpr double GetDistanceFromSurfaceInMeters(void) const noexcept
        {
            return(m_DistanceFromSurfaceInMeters);
        }

        [[nodiscard]] inline constexpr double GetAltitudeInMeters(void) const
        {
            return(GetDistanceFromSurfaceInMeters());
        }

        inline constexpr void Set(double const up_down_angle, double const left_right_angle, double const distance_from_surface) noexcept
        {
            m_UpDownAngleInDegrees = up_down_angle;
            m_LeftRightAngleInDegrees = left_right_angle;
            m_DistanceFromSurfaceInMeters = distance_from_surface;
        }

        inline constexpr void SetUpDownAngleInDegrees(double const up_down_angle) noexcept
        {
            m_UpDownAngleInDegrees = up_down_angle;
        }

        inline constexpr void SetLatitude(double const up_down_angle) noexcept
        {
            SetUpDownAngleInDegrees(up_down_angle);
        }

        inline constexpr void SetLeftRightAngleInDegrees(double const left_right_angle) noexcept
        {
            m_LeftRightAngleInDegrees = left_right_angle;
        }

        inline constexpr void SetLongitude(double const left_right_angle) noexcept
        {
            SetLeftRightAngleInDegrees(left_right_angle);
        }

        inline constexpr void SetDistanceFromSurfaceInMeters(double const distance_from_surface) noexcept
        {
            m_DistanceFromSurfaceInMeters = distance_from_surface;
        }

        inline constexpr void SetAltitudeInMeters(double const distance_from_surface) noexcept
        {
            SetDistanceFromSurfaceInMeters(distance_from_surface);
        }

        inline CPolarCoordinate& operator = (CPolarCoordinate const& source)
        {
            Copy(source);
            return(*this);
        }
    };

    class CEarth
    {
        public:

        enum class Ellipsoid
        {
            Unknown,
            Perfect_Sphere,
            Airy,
            Austrailian_National,
            Bessell_1841,
            Bessel_1841_Nambia,
            Clarke_1866,
            Clarke_1880,
            Everest,
            Fischer_1960_Mercury,
            Fischer_1968,
            GRS_1967,
            GRS_1980,
            Helmert_1906,
            Hough,
            International,
            Krassovsky,
            Modified_Airy,
            Modified_Everest,
            Modified_Fischer_1960,
            South_American_1969,
            Topex_Poseidon_Pathfinder_ITRF,
            WGS_60,
            WGS_66,
            WGS_72,
            WGS_84,
            Custom
        };

        enum class Datum
        {
            NAD_27 = (int)Ellipsoid::Clarke_1866,
            Tokyo = (int)Ellipsoid::Bessell_1841,
        };

        protected:

        // These are the magic numbers. They are the "real" data. They are facts,
        // measurements. Everything else about an ellipse is derived (or computed from)
        // these two data items.

        double m_PolarRadiusInMeters{ 0.0 };
        double m_EquatorialRadiusInMeters{ 0.0 };

        // Here be the things that can be derived from our hard data.
        // We compute these values using the two pieces of data that we
        // know about the ellipse.

        double m_Flattening{ 0.0 };
        double m_EccentricitySquared{ 0.0 };

        // Here's stuff specific to the C++ class

        Ellipsoid m_EllipsoidID{ Ellipsoid::Unknown };

        void m_ComputeFlattening(void) noexcept;
        void m_ComputeEccentricitySquared(void) noexcept;
        void m_Initialize(void) noexcept;

        public:

        CEarth(Ellipsoid const ellipsoid = Ellipsoid::WGS_84);
        CEarth(CEarth const& source);
        ~CEarth() = default;

        void AddLineOfSightDistanceAndDirectionToCoordinate(CPolarCoordinate const& point_1, double const distance, double const direction, CPolarCoordinate& point_2, double const height_above_surface_of_point_2 = 0.0) noexcept;
        void AddSurfaceDistanceAndDirectionToCoordinate(CEarthCoordinate const& point_1, double const distance, double const direction, CEarthCoordinate& point_2) noexcept;
        void AddSurfaceDistanceAndDirectionToCoordinate(CEarthCoordinate const& point_1, double const distance, double const direction, CPolarCoordinate& point_2) noexcept;
        void AddSurfaceDistanceAndDirectionToCoordinate(CPolarCoordinate const& point_1, double const distance, double const direction, CEarthCoordinate& point_2) noexcept;
        void AddSurfaceDistanceAndDirectionToCoordinate(CPolarCoordinate const& point_1, double const distance, double const direction, CPolarCoordinate& point_2) noexcept;
        void Convert(CEarthCoordinate const& cartesian_coordinate, CPolarCoordinate& polar_coordinate) const noexcept;
        void Convert(CPolarCoordinate const& polar_coordinate, CEarthCoordinate& cartesian_coordinate) const noexcept;
        [[nodiscard]] inline double GetDistanceToHorizon(CEarthCoordinate const& point_1) const noexcept
        {
            CPolarCoordinate polar_coordinate;

            Convert(point_1, polar_coordinate);

            return(GetDistanceToHorizon(polar_coordinate));
        }

        [[nodiscard]] inline double GetDistanceToHorizon(CPolarCoordinate const& point_1) const noexcept
        {
            // d = ::sqrt( feet ) *  1.144 for nmi
            // optical horizon is 1.317 * sqrt( h );
            // d= ::sqrt( 17 * height_in_meters ); d is in meters

            auto const distance_to_horizon{ GeodesyFoundationClasses::CMath::SquareRoot(17.0 * point_1.GetDistanceFromSurfaceInMeters()) };

            return(distance_to_horizon);
        }

        [[nodiscard]] inline constexpr double GetEquatorialRadiusInMeters(void) const noexcept
        {
            return(m_EquatorialRadiusInMeters);
        }

        [[nodiscard]] inline constexpr double GetPolarRadiusInMeters(void) const noexcept
        {
            return(m_PolarRadiusInMeters);
        }

        [[nodiscard]] double GetLineOfSightDistanceFromCourse(CEarthCoordinate const& current_location, CEarthCoordinate const& point_a, CEarthCoordinate const& point_b) const noexcept;
        [[nodiscard]] double GetLineOfSightDistance(CEarthCoordinate const& point_1, CEarthCoordinate const& point_2) const noexcept;
        [[nodiscard]] double GetLineOfSightDistance(CPolarCoordinate const& point_1, CEarthCoordinate const& point_2) const noexcept;
        [[nodiscard]] double GetLineOfSightDistance(CEarthCoordinate const& point_1, CPolarCoordinate const& point_2) const noexcept;
        [[nodiscard]] double GetLineOfSightDistance(CPolarCoordinate const& point_1, CPolarCoordinate const& point_2) const noexcept;
        [[nodiscard]] double GetSurfaceDistance(CEarthCoordinate const& point_1, CEarthCoordinate const& point_2, double* heading_from_point_1_to_point_2 = nullptr, double* heading_from_point_2_to_point_1 = nullptr) const noexcept;
        [[nodiscard]] double GetSurfaceDistance(CEarthCoordinate const& point_1, CPolarCoordinate const& point_2, double* heading_from_point_1_to_point_2 = nullptr, double* heading_from_point_2_to_point_1 = nullptr) const noexcept;
        [[nodiscard]] double GetSurfaceDistance(CPolarCoordinate const& point_1, CEarthCoordinate const& point_2, double* heading_from_point_1_to_point_2 = nullptr, double* heading_from_point_2_to_point_1 = nullptr) const noexcept;
        [[nodiscard]] double GetSurfaceDistance(CPolarCoordinate const& point_1, CPolarCoordinate const& point_2, double* heading_from_point_1_to_point_2 = nullptr, double* heading_from_point_2_to_point_1 = nullptr) const noexcept;
        void SetEllipsoid(Ellipsoid const ellipsoid) noexcept;
        void SetEllipsoidByRadii(double const equatorial_radius, double const polar_radius) noexcept;
        void SetEllipsoidByEquatorialRadiusAndFlattening(double const equatorial_radius, double const flattening) noexcept;

        void Copy(CEarth const& source) noexcept;
        CEarth& operator = (CEarth const& source);
    };
}
