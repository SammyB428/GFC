/*
Author: Samuel R. Blackburn
Internet: wfc@pobox.com

"You can get credit for something or get it done, but not both."
Dr. Richard Garwin

The MIT License (MIT)

Copyright (c) 1997-2015 Sam Blackburn

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

inline double CMath::AbsoluteValue( const double& value )
{
   return( ::fabs( value ) );
}

inline double CMath::ArcCosine( const double& value )
{
   return( ::acos( value ) );
}

inline double CMath::ArcSine( const double& value )
{
   return( ::asin( value ) );
}

inline double CMath::ArcTangent( const double& value )
{
   return( ::atan( value ) );
}

inline double CMath::ArcTangentOfYOverX( const double& y, const double& x )
{
   return( ::atan2( y, x ) );
}

inline double CMath::Ceiling( const double& value )
{
   return( ::ceil( value ) );
}

inline void CMath::ConvertDecimalDegreesToDegreesMinutesSeconds( double decimal_degrees, double& degrees, double& minutes, double& seconds )
{
   double fractional_part = 0.0;

   double integer_part = 0;

   fractional_part = ::modf( decimal_degrees, &integer_part );

   degrees = integer_part;

   if ( decimal_degrees < 0.0 )
   {
      fractional_part *= (-1.0);
   }

   minutes = fractional_part * 60.0;

   fractional_part = ::modf( minutes, &integer_part );

   minutes = integer_part;

   seconds = fractional_part * 60.0;
}

inline double CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees( double degrees, double minutes, double seconds )
{
   double decimal_degrees = 0.0;

   decimal_degrees = degrees;

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

inline double CMath::ConvertDegreesToRadians( const double& degrees )
{
   double radians           = 0.0;
   double pi_divided_by_180 = CMath::Pi() / 180.0;
   
   radians = degrees * pi_divided_by_180;

   return( radians );
}

inline double CMath::ConvertRadiansToDegrees( const double& radians )
{
   double degrees = 0.0;

   double conversion_factor = 180.0 / CMath::Pi();

   degrees = radians * conversion_factor;

   return( degrees );
}

inline double CMath::Cosine( const double& value )
{
   return( ::cos( value ) );
}

inline double CMath::HyperbolicCosine( const double& value )
{
   return( ::cosh( value ) );
}

inline double CMath::Pi( void )
{
   return( 3.1415926535897932384626433832795028841971693993751058209749445923078164 );
}

inline double CMath::Sine( const double& value )
{
   return( ::sin( value ) );
}

inline double CMath::SquareRoot( const double& value )
{
   return( ::sqrt( value ) );
}