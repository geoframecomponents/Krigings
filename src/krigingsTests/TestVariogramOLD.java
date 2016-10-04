package krigingsTests;


import java.util.HashMap;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import experimentalVariogram.OmsVariogramOLD;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

public class TestVariogramOLD extends HMTestCase {
    @SuppressWarnings("nls")
    public void testVariogram() throws Exception {

        //
        String stationIdField = "ID_P";


		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "/Users/marialaura/Desktop/testW/final_station.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file ="/Users/marialaura/Desktop/testW/P1994_1995.csv";
		reader.idfield = "ID_P";
		reader.tStart = "1994-01-01 00:00";
		reader.tTimestep = 60;
		 reader.tEnd = "1994-01-01 02:00";
		reader.fileNovalue = "-9999";

        reader.initProcess();

        OmsVariogramOLD Meuse = new OmsVariogramOLD();
        Meuse.pm = pm;

        Meuse.inStations = stationsFC;
        Meuse.fStationsid = stationIdField;
        double[][] matrice;

        while( reader.doProcess ) {
            reader.nextRecord();
            HashMap<Integer, double[]> id2ValueMap = reader.outData;
            Meuse.inData = id2ValueMap;
            // Meuse.doPrintfile=true;
            // Meuse.pPath="/Users/giuseppeformetta/Desktop/miofitvariogram/Result.txt";
            Meuse.process();
            /*
             * Extract the result.
             */
            matrice = Meuse.outResult;
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

            }

        }

    }

}
