package VGM;

public class Linear implements Model{
	
	double [] dist;
	double sill;
	double range;
	double nug;
	
	
	public Linear (double[] dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double[] result() {
        int length = dist.length;
        double[] result = new double[length];
        for (int i = 0; i < length; i++) {
            if (dist[i] != 0.0) {
                result[i] = nug + sill * (dist[i] / range);
            }
            if (dist[i] >= range) {
                result[i] = sill;
            }
            //System.out.println(func[i]);
        }
        return result;
	}

}
