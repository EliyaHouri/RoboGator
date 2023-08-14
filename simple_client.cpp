#include <libplayerc++/playerc++.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <unordered_map>
#include <queue>
#include <unordered_set>
#include <map>
#include <thread>
#include <chrono>
#include <algorithm>

using namespace std;

const string VERTEX_FILE_PATH = "csv_files/vertexes.csv";
const string EDGE_FILE_PATH = "csv_files/edges_with_weights_directions.csv";
const string START_POINT_FILE_PATH = "starting_point.txt";

// Global variables for the initial positions in both graph and simulation spaces
double initialSimX = -13.5;
double initialSimY = -4.5;
int initialVertexX = 171;
int initialVertexY = 1256;

// Conversion factors
const double X_factor = 2.925 / 144.0;
const double Y_factor = -5.1 / 205.0;


namespace std {
    template <>
    struct hash<pair<int, int>> {
        size_t operator()(const pair<int, int>& p) const {
            return hash<int>()(p.first) ^ hash<int>()(p.second);
        }
    };
}


// Graph structures
struct Vertex;
struct Edge;

struct Vertex {
    int x, y;
    vector<Edge*> edges;
};

struct Edge {
    Vertex* start;
    Vertex* end;
    int weight;
    string direction;
};

class Graph {
public:
    // Map of vertices with key as (x,y) pair and value as the Vertex object
    unordered_map<pair<int, int>, Vertex*> vertices;
    
    // Add a vertex to the graph
    Vertex* addVertex(int x, int y) {
        Vertex* v = new Vertex{x, y};
        vertices[{x, y}] = v;
        return v;
    }
    
    // Add an edge to the graph
    void addEdge(Vertex* start, Vertex* end, int weight, const string& direction) {
        Edge* e = new Edge{start, end, weight, direction};
        start->edges.push_back(e);
    }
    
    // Get a vertex from the graph given its (x,y) coordinates
    Vertex* getVertex(int x, int y) {
        for (auto& pair : vertices) {
            Vertex* v = pair.second;
            if (v->x == x && v->y == y) {
                return v;
            }
        }
        return nullptr;  // return null if vertex not found
    }


    
    // Dijkstra's algorithm to compute the shortest path between two vertices
    vector<Edge*> shortestPath(Vertex* start, Vertex* end) {
        unordered_set<Vertex*> visited;
        map<Vertex*, Vertex*> previous;
        map<Vertex*, int> distances;

        struct CompareDist {
            map<Vertex*, int>& distancesRef;
            CompareDist(map<Vertex*, int>& distances): distancesRef(distances) {}

            bool operator()(Vertex* a, Vertex* b) {
                return distancesRef[a] > distancesRef[b];
            }
        };

        priority_queue<Vertex*, vector<Vertex*>, CompareDist> pq(distances);

        for (auto& v : this->vertices) {
            distances[v.second] = INT_MAX;  // Initialize distances to "infinity"
            previous[v.second] = nullptr;   // Previous vertex in optimal path from the source
        }

        distances[start] = 0;
        pq.push(start);

        while (!pq.empty()) {
            Vertex* current = pq.top();
            pq.pop();

            if (current == end) {
                break;  // We found the shortest path to the destination
            }

            if (visited.find(current) != visited.end()) {
                continue;  // Skip vertices we've already visited
            }

            visited.insert(current);

            for (Edge* e : current->edges) {
                Vertex* neighbor = (e->start == current) ? e->end : e->start;  // Get the other vertex of the edge
                int newDist = distances[current] + e->weight;
                if (newDist < distances[neighbor]) {
                    distances[neighbor] = newDist;
                    previous[neighbor] = current;
                    pq.push(neighbor);
                }
            }
        }

        // Reconstruct the shortest path
        vector<Edge*> path;
        Vertex* at = end;
        while (at != start) {
            Vertex* prevVertex = previous[at];
            if (!prevVertex) {
                throw runtime_error("Path reconstruction failed!");
            }

            // Find the edge between 'prevVertex' and 'at'
            for (Edge* e : prevVertex->edges) {
                if (e->end == at) {
                    path.insert(path.begin(), e);
                    break;
                }
            }
            at = prevVertex;
        }


        return path;
    }
};


// Global Graph object
Graph graph;



// Helper function to split a string by a delimiter
vector<string> split(const string &s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Helper function to extract x and y from the format (x, y)
pair<int, int> extractCoordinates(const string& token) {
    int x, y;
    sscanf(token.c_str(), "(%d, %d)", &x, &y);
    return {x, y};
}

void loadGraphData() {
    // Load vertices
    {
        int verteces_count = 0;
        ifstream vertex_file("../csv_files/vertexes.csv");
        string line;
        // Skip the header
        getline(vertex_file, line);
        while (getline(vertex_file, line)) {
            vector<string> tokens = split(line, ',');
            if (tokens.size() < 2) {
                continue; // Skip rows with fewer than 2 tokens
            }

            int y = stoi(tokens[0]); // This is the 'Row' value from the CSV
            for (size_t i = 1; i < tokens.size(); ++i) { // Starting from index 1 since index 0 is 'Row' value
                if (tokens[i].empty()) {
                    break; // Move to the next row if the cell is empty
                }

                try {
                    int x = stoi(tokens[i]);
                    graph.addVertex(x, y);
                    verteces_count++;
                } catch(const invalid_argument& e) {
                    // Handle conversion errors gracefully
                    cerr << "Error parsing vertex: x=" << tokens[i] << ", y=" << y << endl;
                }
            }
        }

        vertex_file.close();
    }

    // Load edges
    {
        ifstream edge_file("../csv_files/edges_with_weights_directions.csv");
        string line;
        int edge_count = 0; // count for successfully loaded edges
        // Skip the header
        getline(edge_file, line);
        while (getline(edge_file, line)) {
            vector<string> tokens = split(line, ',');
            try {
                // Remove unwanted characters
                for (size_t i = 0; i < 4; ++i) {
                    tokens[i].erase(remove(tokens[i].begin(), tokens[i].end(), '('), tokens[i].end());
                    tokens[i].erase(remove(tokens[i].begin(), tokens[i].end(), ')'), tokens[i].end());
                    tokens[i].erase(remove(tokens[i].begin(), tokens[i].end(), '"'), tokens[i].end());
                }

                int start_x = stoi(tokens[0]);
                int start_y = stoi(tokens[1]);
                int end_x = stoi(tokens[2]);
                int end_y = stoi(tokens[3]);
                int weight = stoi(tokens[4]);
                string direction = tokens[5];

                Vertex* start = graph.getVertex(start_x, start_y);
                Vertex* end = graph.getVertex(end_x, end_y);
                if (!start || !end) {
                    cerr << "Vertices not found for line: " << line << endl;
                    continue;
                }

                graph.addEdge(start, end, weight, direction);
                edge_count++; // increment the successfully loaded edge count
            } catch(const invalid_argument& e) {
                // Handle conversion errors gracefully
                cerr << "Error parsing line: " << line << ". Error: " << e.what() << endl;
            }
        }
        edge_file.close();
    }
}

double calculateDistance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

bool isRobotStuck(PlayerCc::Position2dProxy& positionProxy, double prevX, double prevY) {
    const double threshold = 0.01;  // A small threshold to determine if the robot has moved
    return calculateDistance(prevX, prevY, positionProxy.GetXPos(), positionProxy.GetYPos()) < threshold;
}

void attemptToFreeRobot(PlayerCc::Position2dProxy& positionProxy) {
    // Turn the robot slightly and move forward a bit
    positionProxy.SetSpeed(0.0, M_PI/6);  // Turn at 30 degrees per second
    sleep(2);  // Turn for 2 seconds
    positionProxy.SetSpeed(0.2, 0.0);  // Move forward
    sleep(2);  // Move for 2 seconds
}

double getMinLaserRange(PlayerCc::LaserProxy& laserProxy) {
    double minRange = std::numeric_limits<double>::max();
    for(int i = 0; i < laserProxy.GetCount(); ++i) {
        if(laserProxy.GetRange(i) < minRange) {
            minRange = laserProxy.GetRange(i);
        }
    }
    return minRange;
}

int avoidObstacle(PlayerCc::Position2dProxy& positionProxy, PlayerCc::LaserProxy& laserProxy, PlayerCc::PlayerClient& robot) {
    const double obstacleThreshold = 0.63;
    const double clearPathThreshold = 0.2;  // Slightly increased for more clearance
    const double rotateSpeed = PlayerCc::dtor(20);  // Rotate at 20 degrees per second
    const double forwardSpeed = 0.2;

    double minRange = getMinLaserRange(laserProxy);

    while (minRange < obstacleThreshold) {
        robot.Read();

        // Determine the best rotation direction (left or right) based on the laser readings.
        // Assuming laserProxy.GetCount() gives the number of readings and laserProxy[i] gives the i-th reading.
        bool rotateLeft = false;
        int halfReadings = laserProxy.GetCount() / 2;
        if (laserProxy[halfReadings] > laserProxy[0]) {
            rotateLeft = true;
        }

        // Rotate in the best direction
        if (!rotateLeft) {
            positionProxy.SetSpeed(0, rotateSpeed);  // Rotate left
        } else {
            positionProxy.SetSpeed(0, -rotateSpeed);  // Rotate right
        }
        sleep(3);

        minRange = getMinLaserRange(laserProxy);

        if (minRange > clearPathThreshold) {
            positionProxy.SetSpeed(forwardSpeed, 0);  // Move forward
            sleep(2.5 * getMinLaserRange(laserProxy));
            positionProxy.SetSpeed(0, 0);  // Stop

            sleep(1);
            return 1;  // Successfully avoided the obstacle
        }
    }
    return 0;  // Obstacle not avoided
}





void navigatePath(PlayerCc::Position2dProxy& positionProxy, vector<Edge*>& path, PlayerCc::PlayerClient& robot, PlayerCc::LaserProxy& laserProxy) {

    const double threshold = 0.05; // Increase the threshold a bit

    for (Edge* edge : path) {
        double targetSimX = initialSimX + (edge->end->x - initialVertexX) * X_factor;
        double targetSimY = initialSimY + (edge->end->y - initialVertexY) * Y_factor;

        while (true) {
            robot.Read();
            // Check for obstacles
            if (getMinLaserRange(laserProxy) < 2) {
                if (avoidObstacle(positionProxy, laserProxy, robot)){
                    positionProxy.GoTo(targetSimX, targetSimY, positionProxy.GetYaw());
                    sleep(2);  // Adjust this value based on how frequently you want to check the position
                }
            }

            double currentSimX = positionProxy.GetXPos();
            double currentSimY = positionProxy.GetYPos();

            double distance = calculateDistance(currentSimX, currentSimY, targetSimX, targetSimY);
            if (distance <= threshold) {
                break;  // The robot is close enough to the target
            } else {
                positionProxy.GoTo(targetSimX, targetSimY, positionProxy.GetYaw());
                sleep(0.2);  // Adjust this value based on how frequently you want to check the position
            }
        }

        // Update the initial positions for the next iteration
        initialSimX = targetSimX;
        initialSimY = targetSimY;
        initialVertexX = edge->end->x;
        initialVertexY = edge->end->y;
    }
}




pair<int, int> getOrSetRandomStartingPoint() {
    return {171, 1256};
}

void logCurrentVertex(const Vertex* v) {
    cout << "Robot is currently at vertex (" << v->x << ", " << v->y << ")." << endl;
}

void logDestVertex(const Vertex* v) {
    cout << "Robot is moving to vertex (" << v->x << ", " << v->y << ")." << endl;
}

void printPath(const vector<Edge*>& path) {
    cout << "Shortest Path:" << endl;
    for (const Edge* edge : path) {
        cout << "(" << edge->start->x << ", " << edge->start->y << ") -> "
             << "(" << edge->end->x << ", " << edge->end->y << ")"
             << " with weight " << edge->weight
             << " and direction " << edge->direction << endl;
    }
}

bool isInputValid(const string& input, int& x, int& y) {
    stringstream ss(input);
    if (!(ss >> x >> y)) {
        return false;  // Format is not correct
    }

    // Check if vertex exists
    Vertex* vertex = graph.getVertex(x, y);
    return vertex != nullptr;
}

void keepConnectionAlive(PlayerCc::PlayerClient& robot) {
    while (true) {
        robot.Read();  // Read from the robot to keep the connection active
        std::this_thread::sleep_for(std::chrono::seconds(1));  // Sleep for 1 second
    }
}

// Function to compute the distance matrix for the given set of points
vector<vector<int>> generateDistanceMatrix(const vector<Vertex*>& inputPoints) {
    int n = inputPoints.size();
    vector<vector<int>> distanceMatrix(n, vector<int>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) continue;  // skip the distance to itself

            // Get the shortest path between the two points
            vector<Edge*> path = graph.shortestPath(inputPoints[i], inputPoints[j]);
            
            // Calculate the distance (sum of weights of the edges in the path)
            int distance = 0;
            for (const auto& edge : path) {
                distance += edge->weight;
            }
            
            distanceMatrix[i][j] = distance;
        }
    }

    return distanceMatrix;
}



// Brute-force TSP solver
vector<int> solveTSP(const vector<vector<int>>& distanceMatrix) {
    int n = distanceMatrix.size();
    vector<int> vertexOrder(n-2);  // Exclude the starting and the last point (meeting point)
    vector<int> bestOrder(n);
    for (int i = 1; i < n-1; ++i) {
        vertexOrder[i-1] = i;
    }

    int shortestPathLength = INT_MAX;

    // Brute-force through all possible permutations
    do {
        int currentPathLength = distanceMatrix[0][vertexOrder[0]];  // starting from predefined starting point
        for (int i = 0; i < n - 3; ++i) {
            currentPathLength += distanceMatrix[vertexOrder[i]][vertexOrder[i + 1]];
        }
        currentPathLength += distanceMatrix[vertexOrder[n-3]][n-1];  // Add distance to the meeting point

        if (currentPathLength < shortestPathLength) {
            shortestPathLength = currentPathLength;
            bestOrder[0] = 0;
            for (int i = 1; i < n-1; ++i) {
                bestOrder[i] = vertexOrder[i-1];
            }
            bestOrder[n-1] = n-1;  // Ensure the last point is the meeting point
        }
    } while (next_permutation(vertexOrder.begin(), vertexOrder.end()));

    return bestOrder;
}


int getMeetingTime() {
    int minutes = 0;
    do {
        cout << "Enter the meeting time in minutes (5-60): ";
        cin >> minutes;
    } while (minutes < 5 || minutes > 60);
    return minutes;
}

Vertex* getMeetingPoint() {
    int x, y;
    cout << "Enter the meeting point, format: x y: ";
    cin >> x >> y;
    while (!isInputValid(to_string(x) + " " + to_string(y), x, y)) {
        cout << "Invalid input or vertex not found. Please try again." << endl;
        cout << "Enter the meeting point, format: x y: ";
        cin >> x >> y;
    }
    return graph.getVertex(x, y);
}

int main(int argc, char *argv[]) {
    using namespace PlayerCc;

    // Initialize random seed
    srand(time(NULL));

    // Initialize player client
    PlayerClient robot("localhost");
    Position2dProxy positionProxy(&robot,0);
    LaserProxy laserProxy(&robot, 0);  // Added LaserProxy


    std::thread connectionThread(keepConnectionAlive, std::ref(robot));
    connectionThread.detach();  // Detach the thread so it runs independently

    loadGraphData();

    // Get the fixed starting point
    auto startingPoint = getOrSetRandomStartingPoint();
    Vertex* current = graph.getVertex(startingPoint.first, startingPoint.second);  // Initialize the current position with the starting point
    logCurrentVertex(current);

    while (true) {
        vector<Vertex*> inputPoints;
        // Add the starting point to the inputPoints
        inputPoints.push_back(current);

        cout << "Enter a list of destination points, format: x y. Type 'done' when finished or 'exit' to quit:" << endl;
        while (true) {
            string input;
            getline(cin, input);

            if (input == "exit") {
                return 0;  // Exit the program
            }

            if (input == "done") {
                if (inputPoints.size() <= 1) {
                    cout << "You didn't enter any points." << endl;
                    continue;  // Continue to prompt the user for more points
                }
                break;  // Exit the loop and proceed to TSP solution
            }


            int x, y;
            stringstream ss(input);
            ss >> x >> y;

            if (!isInputValid(input, x, y)) {
                cout << "Invalid input or vertex not found. Please try again." << endl;
            } else {
                Vertex* destination = graph.getVertex(x, y);
                inputPoints.push_back(destination);
            }
        }
        

        // Get the meeting point and time
        Vertex* meetingPoint = getMeetingPoint();
        inputPoints.push_back(meetingPoint);  // Add meeting point to TSP problem

        int meetingTimeInMinutes = getMeetingTime();

        vector<vector<int>> distanceMatrix = generateDistanceMatrix(inputPoints);
        vector<int> tspOrder = solveTSP(distanceMatrix);

        // Navigate through the TSP path
        cout << "Starting TSP navigation..." << endl;

        int startTime = time(NULL); // Capture the start time

        vector<Edge*> entirePath;  // To store the entire path

        for (int i = 0; i < tspOrder.size() - 1; ++i) {
            Vertex* next = inputPoints[tspOrder[i+1]];
            vector<Edge*> pathSegment = graph.shortestPath(current, next);
            
            // Append the path segment to the entire path
            entirePath.insert(entirePath.end(), pathSegment.begin(), pathSegment.end());

            navigatePath(positionProxy, pathSegment, robot, laserProxy);

            int elapsedTime = time(NULL) - startTime; // Time taken so far
            int timeLeftInSeconds = meetingTimeInMinutes * 60 - elapsedTime;

            if (next == meetingPoint) { // Check if this is the meeting point
                printf("The meeting will start in %02d:%02d:%02d minutes, WELCOME!\n", 
                       timeLeftInSeconds / 3600, (timeLeftInSeconds % 3600) / 60, timeLeftInSeconds % 60);
                break;  // Exit the loop since we've reached the meeting point
            } else if (timeLeftInSeconds > 0) {
                printf("There is %02d:%02d:%02d minutes left for the meeting in (%d, %d).\n", 
                       timeLeftInSeconds / 3600, (timeLeftInSeconds % 3600) / 60, timeLeftInSeconds % 60, 
                       meetingPoint->x, meetingPoint->y);
            }

            // Update current position for next iteration
            current = next;
        }

        // Printing the entire path
        cout << "Hey Gal, Thought you will find interest in the whole path :)" << endl;
        cout << "Entire TSP Path:" << endl;
        printPath(entirePath);
        cout << "------------------------" << endl;

        cout << "Navigation completed!" << endl;
        cout << "------------------------" << endl;

        // Update the starting point for the next iteration to the last visited point (meeting point)
        current = meetingPoint;

        cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // Add this line to clear the buffer at the end of the loop
    }

    return 0;
}
