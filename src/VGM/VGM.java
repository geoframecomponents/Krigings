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
package VGM;


import oms3.annotations.*;

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
public class VGM extends JGTModel {

	@Description("Distances vector.")
	@In
	public double distance;

	@Description("Sill value.")
	@In
	public double sill;

	@Description("Range value.")
	@In
	public double range;

	@Description("Nugget value.")
	@In
	public double nugget;

	@Description("Model name")
	@In
	public String modelname;

	@Description("Semivariance vector.")
	@Out
	public double result;

	Model modelVGM;

	@Execute
	public void process() throws Exception {
		result = calculateVGM(modelname, distance, sill, range, nugget);
	}


	public double calculateVGM( String model,double distance, double sill, double range, double nug) {

		modelVGM=SimpleModelFactory.createModel(model,distance, sill, range, nug);
		double result=modelVGM.result();

		return result;
	}

}
