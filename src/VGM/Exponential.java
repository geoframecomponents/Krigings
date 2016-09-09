package VGM;


public class Exponential implements Model{

	double dist;
	double sill;
	double range;
	double nug;


	public Exponential (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}



	@Override
	public double result() {

		double result = 0;

		if (dist != 0.0) {
			result = nug + sill * (1 - (Math.exp(-dist / range)));
		}
		//System.out.println(func[i]);

		return result;
	}


}
