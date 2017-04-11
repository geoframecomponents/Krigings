/*
 * GNU GPL v3 License
 *
 * Copyright 2016 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
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
package krigingsTests;


import java.util.HashMap;

import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;

import org.junit.Test;
import org.junit.Assert;

import theoreticalVariogram.TheoreticalVariogram;


public class TestTheoreticalVariogram {
	@SuppressWarnings("nls")
	
	@Test
	public void testTheoreticalVariogram() throws Exception {

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file ="resources/Output/experimentalVGM/experimental_distances.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		
		
		OmsTimeSeriesIteratorReader readerV = new OmsTimeSeriesIteratorReader();
		readerV.file ="resources/Output/experimentalVGM/experimental_variogram.csv";
		readerV.idfield = "ID";
		readerV.tStart = "2000-01-01 00:00";
		readerV.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		readerV.fileNovalue = "-9999";

		reader.initProcess();
		readerV.initProcess();

		TheoreticalVariogram Meuse = new TheoreticalVariogram();
		//Meuse.pm = pm;

		Meuse.modelName="power";
        Meuse.range = 2.221328110091;
        Meuse.nugget = 5.050696956032775;
        Meuse.sill= 4.988372025611679;
            
        


		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/theoreticalVGM/theoretical_VGM.csv";
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;


		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			Meuse.inDistanceValues = id2ValueMap;
			
			readerV.nextRecord();
			id2ValueMap = readerV.outData;
			Meuse.inExperimentalVariogramValues= id2ValueMap;


			Meuse.process();

			HashMap<Integer, double[]> resultS = Meuse.outHMtheoreticalVariogram;


			writer.inData = resultS;
			writer.writeNextLine();

		}
		//
		reader.close();
		readerV.close();
		writer.close();



		/*
            for( int i = 0; i < matrice.length; i++ ) {
                double[] expected = new double[3];
                if (i == 0) {
                    expected = new double[]{342, 0.05811439, 450.5313};
                }
                if (i == 1) {
                    expected = new double[]{461, 0.23422428, 825.6636};
                }
                if (i == 2) {
                    expected = new double[]{831, 0.37321886, 716.8723};
                }
                if (i == 3) {
                    expected = new double[]{931, 0.51183712, 742.0123};
                }
                if (i == 4) {
                    expected = new double[]{1022, 0.67333526, 920.7100};
                }
                if (i == 5) {
                    expected = new double[]{1284, 0.81491201, 770.2608};
                }
                if (i == 6) {
                    expected = new double[]{1251, 0.97279054, 757.0037};
                }
                if (i == 7) {
                    expected = new double[]{1663, 1.10785333, 865.5750};
                }
                if (i == 8) {
                    expected = new double[]{1582, 1.26463821, 814.5208};
                }
                if (i == 9) {
                    expected = new double[]{1807, 1.40859854, 851.0686};
                }
                if (i == 10) {
                    expected = new double[]{1738, 1.55417387, 852.2205};
                }
                if (i == 11) {
                    expected = new double[]{1793, 1.71473991, 884.7430};
                }
                if (i == 12) {
                    expected = new double[]{1664, 1.85265247, 1042.8704};
                }
                if (i == 13) {
                    expected = new double[]{1560, 2.00985981, 1030.4563};
                }
                if (i == 14) {
                    expected = new double[]{1490, 2.14546221, 871.4142};
                }
                assertEquals(expected[0], matrice[i][0], 1);
                assertEquals(expected[1], matrice[i][1], 0.1);
                assertEquals(expected[2], matrice[i][2], 1);
		 */

	}

}




