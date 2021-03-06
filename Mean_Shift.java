//http://rsbweb.nih.gov/ij/plugins/download/Mean_Shift.java
import ij.*;
import ij.gui.*;
import ij.plugin.filter.*;
import ij.process.*;
import java.awt.*;

public class Mean_Shift implements ExtendedPlugInFilter, DialogListener {
	int flags = DOES_8G|DOES_16|DOES_32|DOES_RGB|KEEP_PREVIEW|PARALLELIZE_STACKS;
	ImagePlus imp; // ImagePlus object
	int rad, rad2;
	float radCol, radCol2;
	int nPasses, pass;
	boolean isRGB;

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		isRGB = imp!=null && imp.getType()==ImagePlus.COLOR_RGB;
		if (!isRGB)
			flags = flags|CONVERT_TO_FLOAT;
		return flags;
	}

	public void run(ImageProcessor ip) {
		pass++;
		if (isRGB)
			filterRGBImage(ip);
		else
			filterGrayscaleImage(ip);
		showProgress( 1.0);
	}

	public void filterRGBImage(ImageProcessor ip) {
		int width = ip.getWidth();
		int height = ip.getHeight();
		int[] pixels = (int[])ip.getPixels();
		float[][] pixelsf = new float[width*height][3];

		for (int i = 0; i < pixelsf.length; i++) {
			int argb = pixels[i];

			int r = (argb >> 16) & 0xff;
			int g = (argb >>  8) & 0xff;
			int b = (argb) & 0xff;

			pixelsf[i][0] = 0.299f  *r + 0.587f *g + 0.114f  *b; // Y
			pixelsf[i][1] = 0.5957f *r - 0.2744f*g - 0.3212f *b; // I
			pixelsf[i][2] = 0.2114f *r - 0.5226f*g + 0.3111f *b; // Q
		}

		float shift = 0;
		int iters = 0;
		for (int y=0; y<height; y++) {
			if (y%20==0) showProgress( y/(double)height);
			for (int x=0; x<width; x++) {

				int xc = x;
				int yc = y;
				int xcOld, ycOld;
				float YcOld, IcOld, QcOld;
				int pos = y*width + x;
				float[] yiq = pixelsf[pos];
				float Yc = yiq[0];
				float Ic = yiq[1];
				float Qc = yiq[2];

				iters = 0;
				do {
					xcOld = xc;
					ycOld = yc;
					YcOld = Yc;
					IcOld = Ic;
					QcOld = Qc;

					float mx = 0;
					float my = 0;
					float mY = 0;
					float mI = 0;
					float mQ = 0;
					int num=0;

					for (int ry=-rad; ry <= rad; ry++) {
						int y2 = yc + ry;
						if (y2 >= 0 && y2 < height) {
							for (int rx=-rad; rx <= rad; rx++) {
								int x2 = xc + rx;
								if (x2 >= 0 && x2 < width) {
									if (ry*ry + rx*rx <= rad2) {
										yiq = pixelsf[y2*width + x2];

										float Y2 = yiq[0];
										float I2 = yiq[1];
										float Q2 = yiq[2];

										float dY = Yc - Y2;
										float dI = Ic - I2;
										float dQ = Qc - Q2;

										if (dY*dY+dI*dI+dQ*dQ <= radCol2) {
											mx += x2;
											my += y2;
											mY += Y2;
											mI += I2;
											mQ += Q2;
											num++;
										}
									}
								}
							}
						}
					}
					float num_ = 1f/num;
					Yc = mY*num_;
					Ic = mI*num_;
					Qc = mQ*num_;
					xc = (int) (mx*num_+0.5);
					yc = (int) (my*num_+0.5);
					int dx = xc-xcOld;
					int dy = yc-ycOld;
					float dY = Yc-YcOld;
					float dI = Ic-IcOld;
					float dQ = Qc-QcOld;

					shift = dx*dx+dy*dy+dY*dY+dI*dI+dQ*dQ;
					iters++;
				}
				while (shift > 3 && iters < 100);

				int r_ = (int)(Yc + 0.9563f*Ic + 0.6210f*Qc);
				int g_ = (int)(Yc - 0.2721f*Ic - 0.6473f*Qc);
				int b_ = (int)(Yc - 1.1070f*Ic + 1.7046f*Qc);

				pixels[pos] = (0xFF<<24)|(r_<<16)|(g_<<8)|b_;
			}

		}
	}

	public void filterGrayscaleImage(ImageProcessor ip) {
		int width = ip.getWidth();
		int height = ip.getHeight();
		float[] pixelsf = (float[])ip.duplicate().getPixels();
		float[] pixelsResult = (float[])ip.getPixels();

		float shift = 0;
		int iters = 0;
		for (int y=0; y<height; y++) {
			if (y%20==0) showProgress( y/(double)height);
			for (int x=0; x<width; x++) {

				int xc = x;	//start pixel
				int yc = y;	//
				int xcOld, ycOld;
				float YcOld;
				int pos = y*width + x;
				float Yc = pixelsf[pos];//grey scale at (x,y)

				iters = 0;
				do {
					//store current center and grey scale @(xc,yc),value Yc
					xcOld = xc;
					ycOld = yc;
					YcOld = Yc;

					float mx = 0;
					float my = 0;
					float mY = 0;
					int num=0;
					//mean shift smoothing
					//rad : radius(bandwidth) of spatial domain,rad2 : rad^2
					//(xc,yc) is the point to move, so it's a center with radius, rad.
					for (int ry=-rad; ry <= rad; ry++) {
						int y2 = yc + ry;//pts's pos in the radius
						if (y2 >= 0 && y2 < height) {//y2 is valid pos
							for (int rx=-rad; rx <= rad; rx++) {
								int x2 = xc + rx;
								if (x2 >= 0 && x2 < width) {//x2 is in valid pos
									if (ry*ry + rx*rx <= rad2) {//only pts in circle/sphere is considered!

										float Y2 = pixelsf[y2*width + x2];//the grey scale of a point in the sphere

										float dY = Yc - Y2;//delta between center's grey scale to nearby a point

										if (dY*dY <= radCol2) {//radius of range domain (color space)
											//grey scale difference is also in the bandwidth
											//sum up the posXY of the points in bandwidth
											//and sum up the grey scale
											mx += x2;
											my += y2;
											mY += Y2;
											num++;	//number inside spatial radius and range radius
										}
									}
								}
							}
						}
					}
					//points inside the spatial radius are calculated with position and grey scale
					float num_ = 1f/num;
					Yc = mY*num_;			//mean of grey scale within the bandwidth b_s,b_r
					xc = (int) (mx*num_+0.5);	//weigth centerx of the bandwith, circle
					yc = (int) (my*num_+0.5);	//weigth centery of the bandwith
					//generate mean shift vector (dx,dy,dY)
					int dx = xc-xcOld;
					int dy = yc-ycOld;
					float dY = Yc-YcOld;
					//distance of the shift vector
					shift = dx*dx+dy*dy+dY*dY;
					iters++;
				}//shift threshold is "3", if shift >3, move until shift <=3 (converge)
				while (shift > 3 && iters < 100);
				//mean shift converge or exceeds iterations, the mode is found!
				pixelsResult[pos] = Yc;//replace the scale with the mode point's value, smoothing
			}
		}
	}

	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
		GenericDialog gd = new GenericDialog("Mean Shift Filter");
		gd.addNumericField("Spatial Radius: ", 3, 0);
		gd.addNumericField("Color Distance: ", 25, 1);
		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled()) return DONE;
		return IJ.setupDialog(imp, flags);
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		rad = (int) gd.getNextNumber();
		rad2 = rad*rad;
		radCol = (float) (gd.getNextNumber() + 1);
		radCol2 = radCol*radCol;

		return true;
	}

	public void setNPasses(int nPasses) {
		this.nPasses = nPasses;
		pass = 0;
	}

	void showProgress(double percent) {
		percent = (double)(pass-1)/nPasses + percent/nPasses;
		IJ.showProgress(percent);
	}

}



