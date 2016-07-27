package krigings;

public class Linear implements Model{
	
	double dist;
	double sill;
	double range;
	double nug;
	
	
	public Linear (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double  result() {

        double result = 0;

            if (dist != 0.0) {
                result = nug + sill * (dist / range);
            }
            if (dist >= range) {
                result = sill+nug;
            }
 
        return result;
	}

}
