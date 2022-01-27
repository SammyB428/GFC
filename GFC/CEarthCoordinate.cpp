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

#if 0
<HTML>
<HEAD><TITLE>GFC - CEarthCoordinate</TITLE></HEAD>
<BODY>
<H1>CEarthCoordinate</H1>
$Revision: 5 $
<HR>
<H2>Description</H2>
This class encapsulates an Earth-Centered-Earth-Fixed coordinate.
This is also known as a cartesian coordinate. It is made up of three
distances all originating at the center of the Earth.
<H2>Constructors</H2>
<DL COMPACT>
<DT><PRE><B>CEarthCoordinate</B>()
<B>CEarthCoordinate</B>( const CEarthCoordinate&amp; source )</PRE><DD>
Constructs an empty coordinate or copies another <B>CEarthCoordinate</B>.
</DL>
<H2>Methods</H2>
<DL COMPACT>
<DT><PRE>void <B>Copy</B>( const CEarthCoordinate&amp; coordinate )</PRE><DD>
Copies the contents of another <B>CEarthCoordinate</B>.
<DT><PRE>void <B>Get</B>( double&amp; x_coordinate, double&amp; y_coordinate, double&amp; z_coordinate )</PRE><DD>
This allows you to retrieve all the data members in one function call.
<DT><PRE>double <B>GetXCoordinateInMeters</B>( void ) const</PRE><DD>
This method returns the X axis coordinate in meters.
Positive values point towards the intersection of the Prime Meridian and the Equator.
<DT><PRE>double <B>GetYCoordinateInMeters</B>( void ) const</PRE><DD>
This method returns the Y axis coordinate in meters.
Positive values point towards the intersection of 90 degrees east of Prime Meridian and the Equator.
<DT><PRE>double <B>GetZCoordinateInMeters</B>( void ) const</PRE><DD>
This method returns the Z axis coordinate in meters.
Positive values point towards the North Pole, negative values towards the South Pole.
<DT><PRE>void <B>Set</B>( double x_coordinate, double y_coordinate, double z_coordinate )</PRE><DD>
This lets you set all of the data members in a single function call.
<DT><PRE>void <B>SetXCoordinateInMeters</B>( double x_coordinate )</PRE><DD>
This method sets the X axis coordinate in meters.
<DT><PRE>void <B>SetYCoordinateInMeters</B>( double y_coordinate )</PRE><DD>
This method sets the Y axis coordinate in meters.
<DT><PRE>void <B>SetZCoordinateInMeters</B>( double z_coordinate )</PRE><DD>
This method sets the Z axis coordinate in meters.
</DL>
<H2>Operators</H2>
<DL COMPACT>
<DT><PRE><B>=</B> ( const CEarthCoordinate&amp; source )</PRE><DD>
Basically calls <B>Copy</B>().
</DL>
<H2>Example</H2>
<PRE><CODE>#include &lt;stdio.h&gt;
#include &lt;GFC.h&gt;
#pragma hdrstop

void main( void )
{
   <A HREF="CPolarCoordinate.htm">CPolarCoordinate</A> here;
   <B>CEarthCoordinate</B> there;

   // here is 39 degrees 12.152 minutes North Latitude, 76 degrees 46.795 minutes West Longitude
   here.SetUpDownAngleInDegrees(     CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees(  39.0, 12.152, 0.0 ) );
   here.SetLeftRightAngleInDegrees(  CMath::ConvertDegreesMinutesSecondsCoordinateToDecimalDegrees( -76.0, 46.795, 0.0 ) );
   here.SetDistanceFromSurfaceInMeters( 1000.0 );

   <A HREF="CEarth.htm">CEarth</A> earth;
   
   earth.Convert( here, there );
}</CODE></PRE>
<HR><I>Copyright, 1998, <A HREF="mailto:wfc@pobox.com">Samuel R. Blackburn</A></I><BR>
$Workfile: CEarthCoordinate.cpp $<BR>
$Modtime: 5/16/00 8:20a $
</BODY>
</HTML>
#endif

