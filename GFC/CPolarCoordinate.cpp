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

#include "GFC.h"
#pragma hdrstop

#if 0
<WFC_DOCUMENTATION>
<HTML>
<HEAD><TITLE>GFC - CPolarCoordinate</TITLE></HEAD>
<BODY>
<H1>CPolarCoordinate</H1>
$Revision: 5 $
<HR>
<H2>Description</H2>
This class encapsulates a polar coordinate. A polar coordinate is
made up of two angles and one distance. The angles originate at the
center of the Earth. One angle measures up and down (latitude) while
the other measures left and right (longitude). The distance is the
altitude above the surface of the Earth. Positive angles approach north
and east while negative values for an angle run south and west.
<H2>Constructors</H2>
<DL COMPACT>
<DT><PRE><B>CPolarCoordinate</B>()
<B>CPolarCoordinate</B>( const CPolarCoordinate&amp; source )</PRE><DD>
Constructs an empty coordinate or copies another <B>CPolarCoordinate</B>.
</DL>
<H2>Methods</H2>
<DL COMPACT>
<DT><PRE>void <B>Copy</B>( const CPolarCoordinate&amp; coordinate )</PRE><DD>
Copies the contents of another <B>CPolarCoordinate</B>.
<DT><PRE>void <B>Get</B>( double&amp; up_down_angle, double&amp; left_right_angle, double&amp; distance )</PRE><DD>
This allows you to retrieve all the data members in one function call.
<DT><PRE>double <B>GetUpDownAngleInDegrees</B>( void ) const</PRE><DD>
This method returns the up/down angle in degrees.
<DT><PRE>double <B>GetLeftRightAngleInDegrees</B>( void ) const</PRE><DD>
This method returns the left/right angle in degrees.
<DT><PRE>double <B>GetDistanceFromSurfaceInMeters</B>( void ) const</PRE><DD>
This method returns the distance from the surface of the ellipsoid.
<DT><PRE>void <B>Set</B>( double up_down_angle, double left_right_angle, double distance )</PRE><DD>
This lets you set all of the data members in a single function call.
<DT><PRE>void <B>SetUpDownAngleInDegrees</B>( double up_down_angle )</PRE><DD>
This method sets the up/down angle.
<DT><PRE>void <B>SetLeftRightAngleInDegrees</B>( double left_right_angle )</PRE><DD>
This method sets the left/right angle.
<DT><PRE>void <B>SetDistanceFromSurfaceInMeters</B>( double distance )</PRE><DD>
This method sets the distance from the surface of the ellipsoid.
</DL>
<H2>Operators</H2>
<DL COMPACT>
<DT><PRE><B>=</B> ( const CPolarCoordinate&amp; source )</PRE><DD>
Basically calls <B>Copy</B>().
</DL>
<H2>Example</H2>
<PRE><CODE>#include &lt;stdio.h&gt;
#include &lt;GFC.h&gt;
#pragma hdrstop

void main( void )
{
   <B>CPolarCoordinate</B> here;
   <B>CPolarCoordinate</B> there;

   // here is 39 degrees 12.152 minutes North Latitude, 76 degrees 46.795 minutes West Longitude
   here.SetUpDownAngleInDegrees(     CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees(  39.0, 12.152, 0.0 ) );
   here.SetLeftRightAngleInDegrees(  CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees( -76.0, 46.795, 0.0 ) );
   here.SetDistanceFromSurfaceInMeters( 1000.0 );

   there = here;
}</CODE></PRE>
<HR><I>Copyright, 1998, <A HREF="mailto:wfc@pobox.com">Samuel R. Blackburn</A></I><BR>
$Workfile: CPolarCoordinate.cpp $<BR>
$Modtime: 5/16/00 8:21a $
</BODY>
</HTML>
</WFC_DOCUMENTATION>
#endif
