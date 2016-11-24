import java.util.*;
import java.lang.*;
import java.lang.reflect.Array;
import java.io.*;
 
 
class AllPairShortestPath
{
	// shortest distance matrix
	private static int dist[][];
	// preceding node matrix
	private static int pred[][];
	// flow across nodes matrix
	private static int flow[][];
	// 1 denotes nodes which are connected via an edge
	private static int linked[][];
	// minimum flow across an edge on the sneaky path between nodes
	private static int minflow[][];
	// maximum flow across an edge on the sneaky path between nodes
	private static int maxflow[][];
	// average flow per edge across the sneaky path between edges 
	private static double avgflow[][];
	// matrix containing the optimal sneaky path
	private static String pathMatrix[][];
	// denotes infinity
    final static int INF = 9999999;
    // matrix length
    private int MatrixLength;
    // array of parsed traffic flow
    private ArrayList<Paths> tflow = new ArrayList<Paths>();
    // starting point
    private static int sneakyStart;
    // ending point
	private static int sneakyEnd;
	// input file target
	private static String fileExtension = new String("N75");
	// stack to lookup optimal paths
	private static Stack pathStack;
	
	
	// Input:  receives a matrix of ints with known distances between connected nodes
	// Output:  The input matrix will now represent the shortest distance from all nodes to all other nodes.
	//			Additionally, a predecessor matrix is created to aid in path reconstruction:  pred
    void floydWarshall(int graph[][])
    { 
        // initialize the dist, flow, and pred before applying floyd warshall
        for (int i = 0; i < MatrixLength; i++)
            for (int j = 0; j < MatrixLength; j++){
            	// copy edges that exist
            	if (graph[i][j] != 0)
            		dist[i][j] = graph[i][j];
            	// first pass : all zero's become infinity
            	// second pass : zero's without an edge between the nodes become infinity
            	else if (i != j && linked[i][j] == 0)
            		dist[i][j] = INF;
            	// second pass only, path exists between nodes, but no traffic used it.  Stays zero.
            	else
            		dist[i][j] = 0;
            	// initialize flow and pred matrices
                flow[i][j] = 0;
        		pred[i][j] = i+1;
            }


 
        // Floyd Warshall with path reconstruction twist
        for (int k = 0; k < MatrixLength; k++)
        {
            for (int i = 0; i < MatrixLength; i++)
            {
                for (int j = 0; j < MatrixLength; j++)
                {
                    // if node k is on the shortest path between i and j, update the value of dist[i][j]
                	// and then update predecessor matrix to reflect the new path
                    if (dist[i][k] + dist[k][j] < dist[i][j]){
                        dist[i][j] = dist[i][k] + dist[k][j];
                        pred[i][j] = pred[k][j];
                    }
                }
            }
        }
    }
 
    // Input : Matrix to output, writer to print to, and a String to header the block
    void printSolution(int matr[][], PrintWriter writer, String header)
    {
        writer.println(header);
        for (int i=0; i<MatrixLength; ++i)
        {
            for (int j=0; j<MatrixLength; ++j)
            {
            	// useful for printing mid processing Matrices for verification purposes, safe to remove for production
                if (matr[i][j]==INF)
                    writer.print("INF ");
                else
                    writer.print(matr[i][j]+"   ");
            }
            writer.println();
        }
    }
    
    // Input : Beginning and end of the desired path
    // Output : A stack of ints containing the properly ordered nodes of the path
    Stack giveMeThePath(int start, int end){
    	// normalize to java array numbering standards
    	int current = end-1;
    	int begin = start -1;
    	Stack pathStack = new Stack();
    	// begin at the ending node
    	pathStack.push(end);
    	// follow the predecessor path until you arrive at the starting node, adding to the stack
    	while (pred[begin][current] != start){
    		current = pred[begin][current]-1;
    		pathStack.push(current+1);
    	}
    	// cap with the staring node
    	pathStack.push(start);	
		return pathStack;
    }
    
    // Input : receives a list of links in a path, and the traffic using the path
    // Weights the flow matrix to represent the traffic flow on the given path
    void makeAdjacencyMatrix(List<Integer> path, int edgeFlow){
    	int start =  path.get(0);
    	int current;
    	for (int i = 1; i < path.size(); i++){
    		current = path.get(i);
    		flow[current-1][start-1] += edgeFlow;
    		start = current;
    	}
    }
    
    // Input/Output : n/a
    // Creates a matrix of traffic flow from the flow in the input file   
	void parseFlowList() {
		List hereIsThePath;
		for (Paths current : tflow) {
			hereIsThePath = giveMeThePath(current.begin + 1, current.end + 1);
			makeAdjacencyMatrix(hereIsThePath, current.flow);
		}
	}

	// Input : String to denote the file location of the input
	// Output : Matrix of edges and weights connecting nodes from the input file
	// Also populates matrix size, and the target sneaky path, as well as saving
	// the flow for later use
	int[][] parseInput(String filename) {
		MatrixLength = 1000;
		File file = new File(filename);
		BufferedReader reader = null;
		int start = 1;
		int end = 2;
		int weight = 0;
		boolean edge;
		String edgeOrFlow = null;
		int graph[][];
		try {
			// open input file
			reader = new BufferedReader(new FileReader(file));
			String text = null;
			int linecount = 0;
			// read parameters header for matrix size and target sneaky path
			while ((text = reader.readLine()) != null) {
				String[] line = text.split(",");
				if (line.length == 3) {
					MatrixLength = Integer.parseInt(line[0].trim());
					sneakyStart = Integer.parseInt(line[1].trim());
					sneakyEnd = Integer.parseInt(line[2].trim());
					// if the parameters were not at the top of the file, close
					// and reopen to return to the beginning
					if (linecount != 0) {
						reader.close();
						reader = new BufferedReader(new FileReader(file));
					}
					break;
				}
				linecount++;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		// Initialize matrix to size
		graph = new int[MatrixLength][MatrixLength];
		// parse remaining input file between edges and flows
		try {
			String text = null;
			// continue until file end
			while ((text = reader.readLine()) != null) {
				if (!text.equals("")) {
					String[] line = text.split(",");
					// avoid lines with incomplete information
					if (line.length == 4) {
						edgeOrFlow = line[0].trim();
						start = Integer.parseInt(line[1].trim()) - 1;
						end = Integer.parseInt(line[2].trim()) - 1;
						weight = Integer.parseInt(line[3].trim());
					}
					// create edge matrix, verify the nodes in question are
					// within the matrix bounds
					if (edgeOrFlow != null && edgeOrFlow.equals("E") && (start < MatrixLength && start >= 0) && (end < MatrixLength && end >= 0)) {
						graph[start][end] = weight;
					}
					// save flows in a list for later use, verify the nodes in
					// question are within the matrix bounds, and greater than zero
					else if (edgeOrFlow != null && edgeOrFlow.equals("F") && (start < MatrixLength && start >= 0) && (end < MatrixLength && end >= 0) && weight > 0) {
						Paths temp = new Paths(start, end, weight);
						tflow.add(temp);
					}
				}
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		// close file
		} finally {
			try {
				if (reader != null) {
					reader.close();
				}
			} catch (IOException e) {
			}
		}
		// initialize matrix where a 1 exists if the pair of nodes are connected by an edge
		linked = new int[graph.length][graph.length];
		for (int i = 0; i < graph.length; i++) {
			for (int j = 0; j < graph.length; j++) {
				if (graph[i][j] != 0)
					linked[i][j] = 1;
			}
		}
		// initialize matrices to proper size
        dist = new int[MatrixLength][MatrixLength];
        pred = new int[MatrixLength][MatrixLength];
        flow = new int[MatrixLength][MatrixLength];
        minflow = new int[MatrixLength][MatrixLength];
    	maxflow = new int[MatrixLength][MatrixLength];
    	avgflow = new double[MatrixLength][MatrixLength];
    	pathMatrix = new String [MatrixLength][MatrixLength];
    	return graph;
    }
    
	// Input : Shortest Distance matrix and Predecessor matrix
	//         Creates weighted matrix containing traffic flow from every node to every node 
	//         in preparation for second pass of Floyd Warshall to find the sneaky path
    void generateOutput(int[][] input, int[][] predessor){
    	double count;
    	// Iterate over input Matrix
    	for (int i = 0; i < input.length; i++){
    		for (int j = 0; j < input.length; j++){
    	    	int max = 0, min = 0;
    			count = 1;
    			Stack pathStack = new Stack();
    			// Appropriately label diagonals
    			if (i == j)
    				pathMatrix[i][j] = Integer.toString(input[i][j]);
    			else{
    				int currentval;
    		    	int current = j;
    		    	int begin = i ;
    		    	int next = j;
    		    	// push ending point to the stack
    		    	pathStack.push(j+1);
    		    	// continue until beginning is reached
    		    	while (pred[begin][current] != (i+1)){
    		    		current = pred[begin][current] - 1;
    		    		currentval = dist[current][next];
    		    		pathStack.push(current+1);
    		    		count++;
    		    		next = current;
    		    		// keep a running min and max for the path
    		    		if (count == 2 || currentval > max)
    		    			max = currentval;   		    	
    		    		if (count == 2 || currentval < min)
    		    			min = currentval;
    		    	}
    		    	// push start to the stack
    		    	pathStack.push(i+1);
		    		current = pred[begin][current] - 1;
		    		currentval = dist[current][next];
		    		// finalize min and max
		    		if (count == 1 || currentval > max)
		    			max = currentval;   		    	
		    		if (count == 1 || currentval < min)
		    			min = currentval;
    			}
    			// loop over stack and create a readable format
    			String pathStackCorrect = "";
    			while(!pathStack.isEmpty()){
    				pathStackCorrect += pathStack.pop();
    				if (!pathStack.isEmpty())
    					pathStackCorrect += ", ";
    			}
    			// save min/max/avg for each path
    			pathMatrix[i][j] = pathStackCorrect;   			
    			avgflow[i][j] = dist[i][j] / count;
    			minflow[i][j] = min;
    			maxflow[i][j] = max;
    			
    		}
    	}
    } 
    
    // Input : Buffer to write to, Matrix of Strings
    //         Prints Matrix to output buffer in a readable format
    void printPath(String matr[][], PrintWriter writer)
    {
			writer.println("Following matrix shows the shortest " + "path between every pair of vertices");
			for (int i = 0; i < MatrixLength; ++i) {
				for (int j = 0; j < MatrixLength; ++j) {
					if (i == j)
						writer.print(String.format("%20s", Integer.toString(i + 1)));
					else
						writer.print(String.format("%20s", matr[i][j]));
				}
				writer.println();
			}
    }
    
    // Input : Buffer to write to, Matrix of doubles
    //         Prints Matrix to output buffer in a readable form
    void printavg(double matr[][], PrintWriter writer)
    {
        writer.println("Following matrix shows the average "+
                         "distances, per edge, between every pair of vertices");
        for (int i=0; i<MatrixLength; ++i)
        {
            for (int j=0; j<MatrixLength; ++j)
            {
            	if (i == j)
            		writer.print(String.format("%20s", Integer.toString((0))));
            	else
            		writer.print(String.format("%20s", Math.round(matr[i][j]*100.0)/100.0));
            }
            writer.println();
        }
    }
    
    // Input : Buffer to print information to, stack of path nodes in int form
    // 		   Prints path to the buffer in readable form
    void printStackToFile(PrintWriter writer, Stack pathStack)
    {
    	writer.print("[");
    		while(pathStack.size() > 0){
    			writer.print(pathStack.pop());
    			if (pathStack.size() != 0)
    				writer.print(", ");
    		}
    		writer.print("]");
			writer.flush();

    }

	public static void main(String[] args) {
		// start time
		long startTime = System.nanoTime();
		
    	try{
    		// open input file and initialize variables
    		PrintWriter writer = new PrintWriter(fileExtension + ".txt", "UTF-8");
    		Stack thePath;
			ArrayList testFlow = new ArrayList();
			AllPairShortestPath a = new AllPairShortestPath();
			// parse input into a matrix
			int graph[][] = a.parseInput("sneakypathinput" + fileExtension + ".txt");
    		writer.println("Calculating sneaky path from " + sneakyStart + " to " + sneakyEnd + ".");
			// find the solution prior to traffic weighting
			a.floydWarshall(graph);
			// add traffic
			a.parseFlowList();
			// find the sneaky path
			a.floydWarshall(flow);
			a.generateOutput(flow, pred);

			// end time
			long endTime = System.nanoTime();
			System.out.println("Computation took " + ((endTime - startTime) / 1000000) + " milliseconds");
			// get the specific sneaky path
			thePath = (a.giveMeThePath(sneakyStart, sneakyEnd));
			
			// print results to the file.
			writer.println();
			writer.println("Here is that sneakypath.");
			a.printStackToFile(writer, thePath);
			writer.println();
			writer.println();
			writer.println("The edge with the lowest number of other cars in the sneaky path is.");
			writer.println(minflow[sneakyStart-1][sneakyEnd-1]);
			writer.println();
			writer.println("The edge with the highest number of other cars in the sneaky path is.");
			writer.println(maxflow[sneakyStart-1][sneakyEnd-1]);
			writer.println();
			writer.println("The average number of other cars on the sneaky path is.");
			writer.println(avgflow[sneakyStart-1][sneakyEnd-1]);
			writer.flush();			writer.close();
    	}
    	catch(Exception x){
    		x.printStackTrace();
    	}
	}
}