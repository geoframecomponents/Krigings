package VGM;

public class Circolar implements Model{
	
	double [] dist;
	double sill;
	double range;
	double nug;
	
	
	public Circolar (double[] dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double[] result() {
		int length = dist.length;
        double hr;
        double[] result = new double[length];
        for (int i = 0; i < length; i++) {
            hr = dist[i] / (range);
            if (dist[i] != 0.0) {
            	result[i] = nug + sill * ((2.0 / Math.PI) * (hr * Math.sqrt(1.0 - hr * hr) + Math.asin(hr)));
            }
            if (dist[i] >= range) {
            	result[i] = sill+nug;
            }
            //System.out.println(func[i]);

        }
        return result;
	}

}
