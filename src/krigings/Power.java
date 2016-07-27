package krigings;

public class Power implements Model{
	
	double dist;
	double sill;
	double range;
	double nug;
	
	
	public Power (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double result() {
        double result = 0;

            if (dist != 0.0) {
                result = nug + sill * (Math.pow(dist, range));
            }

        return result;
	}

}
