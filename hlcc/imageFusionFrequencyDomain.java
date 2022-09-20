package edu.stanford.rsl.tutorial.hlcc;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid2DComplex;

public class imageFusionFrequencyDomain {
	public float c = 0.2f; //cut-off frequency
	public int startAngle = 0; // in degree
	public int endAngle = 160; // in degree
	public imageFusionFrequencyDomain(){
	
	}
	
	/**
	 * get the frequency components after image fusion
	 * @param reconLimited
	 * @param reconHLCC
	 * @return
	 */
	public Grid2DComplex frequencyAfterImageFusion(Grid2D reconLimited, Grid2D reconHLCC, boolean useCircle){
		Grid2DComplex complexLimited = new Grid2DComplex(reconLimited,false);
		complexLimited.transformForward();
		complexLimited.fftshift();
		Grid2D mag = complexLimited.getMagnSubGrid(0, 0, reconLimited.getSize()[0], reconLimited.getSize()[1]);
		mag.getGridOperator().log(mag);
		mag.show("Fourier Limited, log magnitude");
		
		Grid2DComplex complexHLCC = new Grid2DComplex(reconHLCC, false);
		complexHLCC.transformForward();
		complexHLCC.fftshift();
		Grid2D mag2 = complexHLCC.getMagnSubGrid(0, 0, reconLimited.getSize()[0], reconLimited.getSize()[1]);
		mag2.getGridOperator().log(mag);		
		mag2.show("Fourier HLCC, log magnitude");
		
		Grid2D smoothMask = this.getSmoothWedgeMask(reconLimited, this.startAngle, this.endAngle, this.c, useCircle);
		
		
		smoothMask.show("mask");
		
		for(int i = 0; i < reconLimited.getSize()[0]; i++)
			for(int j = 0; j < reconLimited.getSize()[1]; j++){
				complexLimited.multiplyAtIndex(i, j, smoothMask.getAtIndex(i, j));
				complexHLCC.multiplyAtIndex(i, j, 1 - smoothMask.getAtIndex(i, j));
				complexLimited.setImagAtIndex(i, j, complexLimited.getImagAtIndex(i, j) + complexHLCC.getImagAtIndex(i, j));
				complexLimited.setRealAtIndex(i, j, complexLimited.getRealAtIndex(i, j)+complexHLCC.getRealAtIndex(i, j));
			}	
		return complexLimited;
	}
	
	public Grid2DComplex frequencyFusion(Grid2DComplex complexLimited, Grid2DComplex complex2, Grid2D mask){
		for(int i = 0; i < mask.getSize()[0]; i++)
			for(int j = 0; j < mask.getSize()[1]; j++){
				complexLimited.multiplyAtIndex(i, j, mask.getAtIndex(i, j));
				complex2.multiplyAtIndex(i, j, 1 - mask.getAtIndex(i, j));
				complex2.setImagAtIndex(i, j, complexLimited.getImagAtIndex(i, j) + complex2.getImagAtIndex(i, j));
				complex2.setRealAtIndex(i, j, complexLimited.getRealAtIndex(i, j) + complex2.getRealAtIndex(i, j));
			}	
		return complex2;
	}
	
	public Grid2DComplex getDoubleWedgeArea(Grid2DComplex complexLimited, Grid2D mask){
		Grid2DComplex complex2 = new Grid2DComplex(complexLimited);
		for(int i = 0; i < mask.getSize()[0]; i++)
			for(int j = 0; j < mask.getSize()[1]; j++){
				complex2.setImagAtIndex(i, j, complex2.getImagAtIndex(i, j) * (1-mask.getAtIndex(i, j)));
				complex2.setRealAtIndex(i, j, complex2.getRealAtIndex(i, j) * (1-mask.getAtIndex(i, j)));
			}	
		return complex2;
	}
	
	/**
	 * the circular area where Fourier coefficients can be well estimated
	 * @param mask: double wedge mask
	 * @param radius: the area where can be well estimated
	 * @return
	 */
	void addCircularArea(Grid2D mask, int radius){
		int Xcenter = mask.getWidth()/2;
		int Ycenter = mask.getHeight()/2;
		for(int i = 0; i < radius; i ++)
			for(int j = 0; j < radius; j++){
				if(i*i+j*j<radius*radius){
					mask.setAtIndex(Xcenter+i, Ycenter+j, 0);
					mask.setAtIndex(Xcenter+i, Ycenter-j, 0);
					mask.setAtIndex(Xcenter-i, Ycenter+j, 0);
					mask.setAtIndex(Xcenter-i, Ycenter-j, 0);
				}
			}
		
	}
	
	
	void addCircularArea2(Grid2D mask, int radius){
		int Xcenter = mask.getWidth()/2;
		int Ycenter = mask.getHeight()/2;
		for(int i = 0; i < radius; i ++)
			for(int j = 0; j < radius; j++){
				if(i*i+j*j<radius*radius){
					mask.setAtIndex(Xcenter+i, Ycenter+j, 1);
					mask.setAtIndex(Xcenter+i, Ycenter-j, 1);
					mask.setAtIndex(Xcenter-i, Ycenter+j, 1);
					mask.setAtIndex(Xcenter-i, Ycenter-j, 1);
				}
			}
		
	}
	/**
	 * the zero double wedge region close to the center need to be zero. Because of the Gaussian smoothing, they are not zero now.
	 * @param smoothMask
	 * @param mask
	 * @param radius
	 */
	void dealCenterRegion(Grid2D smoothMask, Grid2D mask, int radius){
		int Xcenter = smoothMask.getWidth()/2;
		int Ycenter = smoothMask.getHeight()/2;
		for(int i = 0; i < radius; i ++)
			for(int j = 0; j < radius; j++){
				if(i*i+j*j<radius*radius){
					smoothMask.setAtIndex(Xcenter+i, Ycenter+j, mask.getAtIndex(Xcenter+i, Ycenter+j));
					smoothMask.setAtIndex(Xcenter+i, Ycenter-j, mask.getAtIndex(Xcenter+i, Ycenter-j));
					smoothMask.setAtIndex(Xcenter-i, Ycenter+j, mask.getAtIndex(Xcenter-i, Ycenter+j));
					smoothMask.setAtIndex(Xcenter-i, Ycenter-j, mask.getAtIndex(Xcenter-i, Ycenter-j));
				}
			}
	}
	
	/**
	 * get the fused image
	 * @param reconLimited
	 * @param reconHLCC
	 * @return
	 */
	public Grid2D getFusedImage(Grid2D reconLimited, Grid2D reconHLCC, boolean useCircle){
		Grid2DComplex complexFused = this.frequencyAfterImageFusion(reconLimited, reconHLCC, useCircle);
		Grid2D mag = complexFused.getMagnSubGrid(0, 0, reconLimited.getSize()[0], reconLimited.getSize()[1]);
	    mag.getGridOperator().log(mag);
		mag.show("frequency magnitude after fusion, logrithmized");
		complexFused.ifftshift();
		complexFused.transformInverse();
		Grid2D imgFused = complexFused.getRealSubGrid(0, 0, reconLimited.getSize()[0], reconLimited.getSize()[1]);
		
		return imgFused;
	}
	
	/**
	 * create a binary double wedge mask for image fusion at the frequency domain
	 * @param img
	 * @param startAngle: angle in degree, e.g. 0
	 * @param endAngle: angle in degree, e.g. 160
	 * @return
	 */
	public Grid2D getWedgeMask(Grid2D img, int startAngle, int endAngle){
		Grid2D mask = new Grid2D(img.getSize()[0], img.getSize()[1]);
		mask.getGridOperator().fill(mask, 1);
		int centerX = img.getSize()[0]/2, centerY = img.getSize()[1]/2;
		float compAngle = 0; //compensation angle
		double slopeStart = Math.tan((startAngle-compAngle)*Math.PI/180.);
		double slopeEnd = Math.tan((endAngle+compAngle)*Math.PI/180.); // is negative
		double ymin, ymax;
		
		for(int x = 0; x <= centerX; x++){
			ymin = (centerX - x)*slopeStart;
			ymax = -(centerX - x)*slopeEnd;
			if(ymin + centerY < 0)
				ymin = - centerY;
			if(ymax + centerY >= img.getSize()[1])
				ymax = centerY-1;
			for(int y = (int) Math.round(ymin); y <= Math.round(ymax); y++){
				mask.setAtIndex(x, y+centerY,0);
			}
		}
		
		for(int x = centerX; x <img.getSize()[0]; x++){
			ymax = (x-centerX)*slopeStart;
			ymin = (x-centerX)*slopeEnd;
			if(ymin + centerY < 0)
				ymin = - centerY;
			if(ymax + centerY >= img.getSize()[1])
				ymax = centerY-1;
			for(int y = (int) Math.round(ymin); y <=Math.round(ymax); y++){
				mask.setAtIndex(x,y+centerY, 0);
			}
		}
		return mask;
	}

	/**
	 * create a smooth double wedge mask to get smooth transition at the wedge boundaries, romove a circur area at center
	 * @param img
	 * @param startAngle
	 * @param endAngle
	 * @param c
	 * @return
	 */
	public Grid2D getSmoothWedgeMask(Grid2D img, int startAngle, int endAngle, float c, boolean useCircle){
		Grid2D mask = this.getWedgeMask(img, startAngle, endAngle);
		if(useCircle)
			this.addCircularArea(mask, 7);
		mask.clone().show("binary mask");
		
		
		Grid2D smoothMask = this.GaussianFilter(mask, c);
		if(useCircle)
			this.addCircularArea(smoothMask, 5);
		//this.dealCenterRegion(smoothMask, mask, 4);
		smoothMask.getGridOperator().removeNegative(smoothMask);
		smoothMask.getGridOperator().setMax(smoothMask, 1);
		//smoothMask.setAtIndex(smoothMask.getWidth()/2, smoothMask.getHeight()/2, 0);//center
		
		return smoothMask;
	}
	
	/**
	 * create a smooth double wedge mask to get smooth transition at the wedge boundaries, at the center add a circular area
	 * @param img
	 * @param startAngle
	 * @param endAngle
	 * @param c
	 * @return
	 */
	public Grid2D getSmoothWedgeMask2(Grid2D img, int startAngle, int endAngle, float c, boolean useCircle){
		Grid2D mask = this.getWedgeMask(img, startAngle, endAngle);
		if(useCircle)
			this.addCircularArea2(mask, 7);
		mask.clone().show("binary mask");
		
		
		Grid2D smoothMask = this.GaussianFilter(mask, c);
		if(useCircle)
			this.addCircularArea2(smoothMask, 5);
		//this.dealCenterRegion(smoothMask, mask, 4);
		smoothMask.getGridOperator().removeNegative(smoothMask);
		smoothMask.getGridOperator().setMax(smoothMask, 1);
		//smoothMask.setAtIndex(smoothMask.getWidth()/2, smoothMask.getHeight()/2, 0);//center
		
		return smoothMask;
	}
	
	/**
	 * Gaussian filter with cut-off frequency c
	 * @param img
	 * @param c
	 * @return
	 */
	public Grid2D GaussianFilter(Grid2D img, float c){
		Grid2DComplex cCopy  = new Grid2DComplex(img,false);		
		cCopy.transformForward();
		//cCopy.fftshift();
		Grid2D kernel = this.GaussianKernel2D(img.getSize(),c);
		for(int i =0; i < img.getSize()[0]; i++)
			for(int j = 0; j < img.getSize()[1]; j++)
				cCopy.multiplyAtIndex(i, j, kernel.getAtIndex(i, j));
		
		//cCopy.ifftshift();
		cCopy.transformInverse();
		Grid2D imgFiltered = cCopy.getRealSubGrid(0, 0, img.getSize()[0], img.getSize()[1]);
		
		return imgFiltered;
	}
	
	/**
	 * compute a 2-D Gaussian kernel with cut-off frequency c, without fftshift
	 * @param size
	 * @param c
	 * @return
	 */
	public Grid2D GaussianKernel2D(int[] size, float c){
		Grid2D kernel = new Grid2D(size[0], size[1]);
		double  x,y,val;
		double d = size[0]*c;
		d = d*d;
		
		for(int i = 0; i< size[0]; i ++){
			for(int j = 0; j < size[1]; j++){
				if(i > size[0]/2)
					x = i -size[0];
				else
					x = i;
				if(j > size[1]/2)
					y = j -size[1];
				else
					y = j;
				
				val = Math.exp(-(x*x+y*y)/(2*d));
				kernel.setAtIndex(i, j,(float) (val));
			}
		}
		
		return kernel;
	}
	
	/**
	 * compute a 2-D Gaussian kernel with cut-off frequency c, with fftshift
	 * @param size
	 * @param c
	 * @return
	 */
	public Grid2D GaussianKernel2DFFTShifted(int[] size, float c){
		Grid2D kernel = new Grid2D(size[0], size[1]);
		double  x,y,val;
		double d = size[0]*c;
		d = d*d;
		
		for(int i = 0; i< size[0]; i ++){
			for(int j = 0; j < size[1]; j++){
					x = i -size[0]/2;
					y = j -size[1]/2;				
				val = Math.exp(-(x*x+y*y)/(2*d));
				kernel.setAtIndex(i, j,(float) (val));
			}
		}
		
		return kernel;
	}
}
