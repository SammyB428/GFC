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

/* SPDX-License-Identifier: MIT */

inline constexpr [[nodiscard]] double GeodesyFoundationClasses::CMath::Pi(void) noexcept
{
    return(3.1415926535897932384626433832795028841971693993751058209749445923078164);
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::AbsoluteValue( double const value ) noexcept
{
   return( std::fabs( value ) );
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::ArcCosine( double const value ) noexcept
{
   return( std::acos( value ) );
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::ArcSine( double const value ) noexcept
{
   return( std::asin( value ) );
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::ArcTangent( double const value ) noexcept
{
   return( std::atan( value ) );
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::ArcTangentOfYOverX( double const y, double const x ) noexcept
{
   return( std::atan2( y, x ) );
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::Ceiling( double const value ) noexcept
{
   return( std::ceil( value ) );
}

inline void GeodesyFoundationClasses::CMath::ConvertDecimalDegreesToDegreesMinutesSeconds( double const decimal_degrees, double& degrees, double& minutes, double& seconds ) noexcept
{
    double integer_part{ 0.0 };

   auto fractional_part = std::modf( decimal_degrees, &integer_part );

   degrees = integer_part;

   if ( decimal_degrees < 0.0 )
   {
      fractional_part *= (-1.0);
   }

   minutes = fractional_part * 60.0;

   fractional_part = std::modf( minutes, &integer_part );

   minutes = integer_part;

   seconds = fractional_part * 60.0;
}

inline constexpr [[nodiscard]] double GeodesyFoundationClasses::CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees( double const degrees, double const minutes, double const seconds ) noexcept
{
   double decimal_degrees{ degrees };

   if ( decimal_degrees < 0.0 )
   {
      decimal_degrees *= (-1.0);
   }

   decimal_degrees += (double) ( minutes / 60.0 );
   decimal_degrees += (double) ( seconds / 3600.0 );

   if ( degrees < 0.0 )
   {
      decimal_degrees *= (-1.0);
   }

   return( decimal_degrees );
}

inline constexpr [[nodiscard]] double GeodesyFoundationClasses::CMath::ConvertDegreesToRadians( double const degrees ) noexcept
{
   double constexpr pi_divided_by_180{ GeodesyFoundationClasses::CMath::Pi() / 180.0 };
   
   auto const radians{ degrees * pi_divided_by_180 };

   return( radians );
}

inline constexpr [[nodiscard]] double GeodesyFoundationClasses::CMath::ConvertRadiansToDegrees( double const radians ) noexcept
{
   auto constexpr conversion_factor{ 180.0 / GeodesyFoundationClasses::CMath::Pi() };

   auto const degrees{ radians * conversion_factor };

   return( degrees );
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::Cosine( double const value ) noexcept
{
   return( ::cos( value ) );
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::HyperbolicCosine( double const value ) noexcept
{
   return( ::cosh( value ) );
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::Sine( double const value ) noexcept
{
   return( std::sin( value ) );
}

inline [[nodiscard]] double GeodesyFoundationClasses::CMath::SquareRoot( double const value ) noexcept
{
   return( std::sqrt( value ) );
}