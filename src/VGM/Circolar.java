package VGM;

public class Circolar implements Model{
	
	double dist;
	double sill;
	double range;
	double nug;
	
	
	public Circolar (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double result() {

        double hr;
        double  result = 0;
            hr = dist / range;
            if (dist!= 0.0) {
            	result = nug + sill * ((2.0 / Math.PI) * (hr * Math.sqrt(1.0 - hr * hr) + Math.asin(hr)));
            } else if (dist >= range) {
            	result  = sill+nug;
            } else if (dist==0){
            	result=0;
            }

        return result;
	}

}
