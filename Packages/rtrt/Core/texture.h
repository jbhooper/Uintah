#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include <GL/glu.h>
#include <GL/glx.h>
#include <math.h>
#define textureHeight 128
#define textureWidth 128

namespace rtrt {
  
template <class T>
class Texture
{
 public:
  Texture()
    {
      //printf( "In Texture creation\n" ); 
      int i, j, k;
      for( i = 0; i < textureHeight; i++ )
	for( j = 0; j < textureWidth; j++ )
	  for( k = 0; k < 4; k++ )
	    textArray[i][j][k] = 0;
      current_color[0] = 0;
      current_color[1] = 1;
      current_color[2] = 0;
      colormap_x_offset = 160;
      colormap_y_offset = 215;
    }

  T textArray[textureHeight][textureWidth][4];
  GLuint TextName;
  float current_color[3];
  int colormap_x_offset;
  int colormap_y_offset;

  /*********************************/
  // COLORMAP CODE
  /*********************************/
  void 
    assign_color(float color[3], float r, float g, float b) 
    {
      //printf( "In Texture::assign_color\n" );
      color[0] = r;
      color[1] = g;
      color[2] = b;
    }


  void 
    interpolate_color(float color_1[3], float color_2[3],
		      float out_color[3], float interpolate) 
    {
      //printf( "In Texture::interpolate_color\n" );
      int i;
      /*  printf("interpolate = %g\n", interpolate);*/
      if (interpolate < 0) interpolate = 0;
      else if (interpolate > 1) interpolate = 1;
  
      for(i = 0; i < 3; i++) 
	{
	  out_color[i] = color_1[i] * (1-interpolate) + color_2[i] * interpolate;
	  /*    printf("i = %d, out_color = %g, color_1 = %g, color_2 = %g\n", i,
		out_color[i], color_1[i], color_2[i]);*/
	}
    }

  void 
    colormap(int x, int y, int dx, int dy, float color[3]) 
    {
      //printf( "In Texture::colormap\n" );
      const int xmin = 0;
      const int ymin = 30;
      const int xmax = 400;
      const int ymax = 400;
      x = (colormap_x_offset + dx + xmax)%xmax;
      y = colormap_y_offset + dy;
      if( y > ymax )
	y = ymax;
      if( y < ymin )
	y = ymin;
      float bottom_hue[3], top_hue[3];
      float hue_index;
      float hue_interpolant;
      float hue_size = (float)1/6;
      int ymid = (ymax - ymin)/2 + ymin;

      /* find the hue */
      hue_index = ((float)x-xmin)/(xmax-xmin);
      /* Since there are only six hues to check we'll do a linear search.  It
	 saves on programming time. */
  
      if (hue_index < hue_size) 
	{
	  /* from red to yello */
	  assign_color(bottom_hue, 1,0,0);
	  assign_color(top_hue, 1,1,0);
	  hue_interpolant = hue_index/hue_size;
	}
      else if (hue_index < hue_size*2)
	{
	  /* from yello to green */
	  assign_color(bottom_hue, 1,1,0);
	  assign_color(top_hue, 0,1,0);
	  hue_interpolant = (hue_index-hue_size)/hue_size;
	}
      else if (hue_index < hue_size*3)
	{
	  /* from green to cyan */
	  assign_color(bottom_hue, 0,1,0);
	  assign_color(top_hue, 0,1,1);
	  hue_interpolant = (hue_index-hue_size*2)/hue_size;
	} 
      else if (hue_index < hue_size*4)
	{
	  /* from cyan to blue*/
	  assign_color(bottom_hue, 0,1,1);
	  assign_color(top_hue, 0,0,1);
	  hue_interpolant = (hue_index-hue_size*3)/hue_size;
	}
      else if (hue_index < hue_size*5) 
	{
	  /* from blue to magenta */
	  assign_color(bottom_hue, 0,0,1);
	  assign_color(top_hue, 1,0,1);
	  hue_interpolant = (hue_index-hue_size*4)/hue_size;
	} 
      else
	{
	  /* from magenta to red*/
	  assign_color(bottom_hue, 1,0,1);
	  assign_color(top_hue, 1,0,0);
	  hue_interpolant = (hue_index-hue_size*5)/hue_size;
	}

      interpolate_color(bottom_hue, top_hue, color, hue_interpolant);

      /* Now to do the interpolation with the black and white components. */

      /*  printf("y = %d, y(%d, %d, %d)\n", y, ymin, ymid, ymax);*/
      if (y < ymid)
	{
	  /* from black to hue */
	  assign_color(bottom_hue, 0,0,0);
	  assign_color(top_hue, color[0], color[1], color[2]);
	  hue_interpolant = ((float)y-ymin)/(ymid-ymin);
	} 
      else
	{
	  /* from hue to white */
	  assign_color(bottom_hue, color[0], color[1], color[2]);
	  assign_color(top_hue, 1,1,1);
	  hue_interpolant = ((float)y-ymid)/(ymax-ymid);
	}

      interpolate_color(bottom_hue, top_hue, color, hue_interpolant);

      for( int i = 0; i < textureHeight; i++ )
	for( int j = 0; j < textureWidth; j++ )
	  {
	    textArray[i][j][0] = color[0];
	    textArray[i][j][1] = color[1];
	    textArray[i][j][2] = color[2];
	  }

      colormap_x_offset = x;
      colormap_y_offset = y;
    }
  /*********************************/
  // END OF COLORMAP CODE
  /*********************************/

  void
    makeOneDimTextureImage( void )
    {
      //printf( "In Texture::makeOneDimTextureImage\n" );
      int i, j;
      float intensity;
      for( i = 0; i < textureHeight; i++ )
	for( j = 0; j < textureWidth; j++ )
	  {
	    intensity = (textureWidth/2-fabs(j-textureWidth/2.0f))/(textureWidth/2);
	    textArray[i][j][0] = current_color[0];
	    textArray[i][j][1] = current_color[1];
	    textArray[i][j][2] = current_color[2];
	    textArray[i][j][3] = intensity;
	  }
			
      glPixelStoref( GL_UNPACK_ALIGNMENT, 1 );
      glGenTextures( 1, &TextName );
      glBindTexture( GL_TEXTURE_2D, TextName );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
      glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGBA, GL_FLOAT, textArray );
    }



  void
    makeEllipseTextureImage( void )
    {
      //printf( "In Texture::makeEllipseTextureImage\n" );
      int i, j;
      float intensity;
      for( i = 0; i < textureHeight; i++ )
	for( j = 0; j < textureWidth; j++ )
	  {
	    intensity = 1.0f - 5*sqrt( (j-textureHeight/2.0f)*(j-textureHeight/2.0f) + (i-textureWidth/2.0f)*(i-textureWidth/2.0f) )/
	                      sqrt( textureHeight*textureHeight/4.0f + textureWidth*textureWidth/4.0f );
	    if( intensity < 0 )
	      intensity = 0;
	    textArray[i][j][0] = current_color[0];
	    textArray[i][j][1] = current_color[1];
	    textArray[i][j][2] = current_color[2];
	    textArray[i][j][3] = intensity;
	  }

      glPixelStoref( GL_UNPACK_ALIGNMENT, 1 );
      glGenTextures( 1, & TextName );
      glBindTexture( GL_TEXTURE_2D, TextName );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
      glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGBA, GL_FLOAT, textArray );
    }



  void
    makeDefaultTextureImage( void )
    {
      //printf( "In Texture::makeDefaultTextureImage\n" );
      int i, j;
      float red = 1;
      float green = 0;
      float blue = 0;
      float hue_width = textureWidth/6;
      float intensity;
      float color_step = 1/hue_width;
      for( i = 0; i < textureHeight; i++ )
	{
	  red = 1;
	  green = blue = 0;
	  intensity = i*100/textureHeight;
	  for( j = 0; j < textureWidth; j++ )
	    {
	      if( j < hue_width )
		green += color_step;
	      else if( j < 2*hue_width )
		red -= color_step;
	      else if( j < 3*hue_width )
		blue += color_step;
	      else if( j < 4*hue_width )
		green -= color_step;
	      else if( j < 5*hue_width )
		red += color_step;
	      else if( j < 6*hue_width )
		blue -= color_step;

	      if( red < 0.0f )
		red = 0.0f;
	      if( green < 0.0f )
		green = 0.0f;
	      if( blue < 0.0f )
		blue = 0.0f;

	      textArray[i][j][0] = red;
	      textArray[i][j][1] = green;
	      textArray[i][j][2] = blue;
	      textArray[i][j][3] = intensity;
	    }
	}			
			
      glPixelStoref( GL_UNPACK_ALIGNMENT, 1 );
      glGenTextures( 1, & TextName );
      glBindTexture( GL_TEXTURE_2D, TextName );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
      glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
      glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGBA, GL_FLOAT, textArray );
    }
};

} // end namespace rtrt

#endif
