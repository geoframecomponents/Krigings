package VGM;

public class Gaussian implements Model{

	double dist;
	double sill;
	double range;
	double nug;


	public Gaussian (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}



	@Override
	public double result() {

		double result=0;
		double hr= dist / (range);

		if (dist != 0.0) {
			 result= nug + sill * (1.0 - (Math.exp(-(hr * hr))));
		}

		return result;
	}


}
