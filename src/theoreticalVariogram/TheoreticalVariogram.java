/* This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package theoreticalVariogram;


import java.util.HashMap;
import java.util.Set;
import java.util.Map.Entry;

import oms3.annotations.*;

import org.geotools.feature.SchemaException;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;




@Description("Teorethical semivariogram models.")
@Documentation("vgm.html")
@Author(name = "Giuseppe Formetta, Adami Francesco & Marialaura Bancheri", contact = " http://www.ing.unitn.it/dica/hp/?user=rigon")
@Keywords("Kriging, Hydrology")
@Label(JGTConstants.STATISTICS)
@Name("kriging")
@Status(Status.CERTIFIED)
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class TheoreticalVariogram extends JGTModel {

	@Description("Distances input Hashmap")
	@In
	public HashMap<Integer, double[]> inDistanceValues;

	@Description("Experimental Variogram input Hashmap")
	@In
	public HashMap<Integer, double[]> inExperimentalVariogramValues;
	
	
	@Description("Distances value.")
	@In
	public boolean doCalibrate;

	@Description("Distances value.")
	@In
	public double distance;

	@Description("Sill value.")
	@In
	@Out
	public double sill;

	@Description("Range value.")
	@In
	@Out
	public double range;

	@Description("Nugget value.")
	@In
	@Out
	public double nugget;

	@Description("Model name")
	@In
	@Out
	public String modelName;

	@Description("the output hashmap withe the semivariance")
	@Out
	public double [] result;

	@Description("the output hashmap withe the semivariance")
	@Out
	public double [] observation;

	@Description("the output hashmap withe the semivariance")
	@Out
	public HashMap<Integer, double[]> outHMtheoreticalVariogram=new HashMap<Integer, double[]>();;

	Model modelVGM;

	@Execute
	public void process() throws Exception {

		// reading the ID of all the stations 
		Set<Entry<Integer, double[]>> entrySet = inDistanceValues.entrySet();
		result=  new double[inDistanceValues.size()];
		observation=  new double[inDistanceValues.size()];

		for (Entry<Integer, double[]> entry : entrySet) {
			Integer ID = entry.getKey();

			distance=inDistanceValues.get(ID)[0];		

			result[ID]=calculateVGM(modelName, distance, sill, range, nugget);
			observation[ID]=inExperimentalVariogramValues.get(ID)[0];
		
			storeVariogramResults(ID,result[ID]);
		}


	}





	public double calculateVGM( String model,double distance, double sill, double range, double nug) {

		modelVGM=SimpleModelFactory.createModel(model,distance, sill, range, nug);
		double result=modelVGM.result();

		return result;
	}

	
	private void storeVariogramResults(int ID, double result) 
			throws SchemaException {

		outHMtheoreticalVariogram.put(ID, new double[]{result});
		
	}
}

