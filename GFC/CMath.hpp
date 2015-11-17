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

#if ! defined( MATH_CLASS_HEADER_FILE )

#define MATH_CLASS_HEADER_FILE

#define HALF_OF_PI          1.5707963267948966
#define ONE_FOURTH_OF_PI    0.78539816339744833
#define PI                  3.14159265358979323846
#define TWO_TIMES_PI        6.2831853071795864769
#define RADIANS_TO_DEGREES 57.29577951308232
#define DEGREES_TO_RADIANS   .0174532925199432958

class CMath
{
   public:

      CMath();
      virtual ~CMath();

      double Cosine(  const double value );
      double Tangent( const double value );
};

#endif // EARTH_CLASS_HEADER_FILE
