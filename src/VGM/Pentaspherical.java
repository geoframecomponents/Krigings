package VGM;

public class Pentaspherical implements Model{
	
	double dist;
	double sill;
	double range;
	double nug;
	
	
	public Pentaspherical (double dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double result() {
    
        double result = 0;
        double hr = 0, h2r2;
            hr = dist / (range);
            h2r2 = hr * hr;
            if (dist != 0.0) {
                result = nug + sill * (hr * ((15.0 / 8.0) + h2r2 * ((-5.0 / 4.0) + h2r2 * (3.0 / 8.0))));
            }
            if (dist >= range) {
                result = sill;
            }
            //System.out.println(func[i]);

        return result;
	}

}
