#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>



const unsigned int number_of_objectives = 1;
const int source_vertex = 0;

struct Edge
{
public:
    int source;
    int destination;
    int* objectives;
    Edge()
    {
        this->source = -1;
        this->destination = -1;
        this->objectives = new int[number_of_objectives];
    }
    Edge(int source, int destination)
    {
        this->source = source;
        this->destination = destination;
        this->objectives = new int[number_of_objectives];
    }
    Edge(const Edge& edge)
    {
        this->source = edge.source;
        this->destination = edge.destination;
        for (unsigned int i = 0; i < number_of_objectives; ++i)
        {
            this->objectives[i] = edge.objectives[i];
        }
    }
    Edge& operator= (const Edge& edge)
    {
        this->source = edge.source;
        this->destination = edge.destination;
        for (unsigned int i = 0; i < number_of_objectives; ++i)
        {
            this->objectives[i] = edge.objectives[i];
        }
        return *this;
    }
    ~Edge()
    {
        this->source = -1;
        this->destination = -1;
        delete[] this->objectives;
        this->objectives = nullptr;
    }
};

struct Couple
{
public:
    int parent;
    int cost;
    Couple()
    {
        this->parent = -1;
        this->cost = -1;
    }
    Couple(int parent, int cost)
    {
        this->parent = parent;
        this->cost = cost;
    }
};

std::string read_graph_from_file(const char* file_name)
{
    std::string graph_file;
    std::ifstream file(file_name);
    if (!file)
    {
        std::cerr << "Failed to open the file " << file_name << ".\n";
        return graph_file;
    }
    std::string line;
    while (std::getline(file, line))
    {
        graph_file += line + '\n';
    }
    file.close();
    return graph_file;
}

void create_edge_list(const std::string& graph_file, std::vector<Edge*>& edge_list)
{
    std::cout << graph_file << '\n';
    std::string number = "";
    unsigned int comma_count = 0;
    unsigned int objective_count = 0;
    Edge* edge = new Edge();
    for (unsigned int i = 0; i < graph_file.length(); ++i)
    {
        if (graph_file[i] >= '0' && graph_file[i] <= '9')
        {
            number += graph_file[i];
        }
        else if (graph_file[i] == ',')
        {
            if (comma_count == 0)
            {
                edge->source = std::stoi(number);
                number = "";
                ++comma_count;
            }
            else if (comma_count == 1)
            {
                edge->destination = std::stoi(number);
                number = "";
                ++comma_count;
            }
            else if (comma_count > 1)
            {
                edge->objectives[objective_count] = std::stoi(number);
                number = "";
                ++objective_count;
            }
            continue;
        }
        else if (graph_file[i] == '\n')
        {
            edge_list.push_back(edge);
            comma_count = 0;
            objective_count = 0;
            edge = new Edge();
            continue;
        }
        else if (graph_file[i] == '~')
        {
            edge_list.push_back(edge);
            comma_count = 0;
            objective_count = 0;
            break;
        }
        else
        {
            std::cerr << "Invalid character found in graph file.\n";
            return;
        }
    }
    return;
}

int find_minimum(std::unordered_map<int, Couple>& vertices)
{
    int minimum = INT_MAX;
    int vertex = -1;
    for (auto& pair : vertices)
    {
        if (pair.second.cost < minimum)
        {
            minimum = pair.second.cost;
            vertex = pair.first;
        }
    }
    return vertex;
}

void make_inital_shortest_path_tree(std::vector<Edge*>& edge_list, std::unordered_map<int, Couple>& shortest_path_tree)
{
    std::unordered_map<int, Couple> not_visited;
    for (Edge* edge : edge_list)
    {
        not_visited[edge->source] = Couple(-1, INT_MAX);
        not_visited[edge->destination] = Couple(-1, INT_MAX);
    }
    not_visited[source_vertex].cost = 0;
    while (!not_visited.empty())
    {
        int current = find_minimum(not_visited);
        if (current == -1)
        {
            break;
        }
        Couple current_data = not_visited[current];
        shortest_path_tree[current] = current_data;
        not_visited.erase(current);
        if (current_data.cost == INT_MAX)
        {
            continue;
        }
        for (unsigned int i = 0; i < edge_list.size(); ++i)
        {
            if (edge_list[i]->source != current)
            {
                continue;
            }
            int neighbor = edge_list[i]->destination;
            if (not_visited.find(neighbor) != not_visited.end())
            {
                int weight = edge_list[i]->objectives[0];
                int new_cost = current_data.cost + weight;
                if (new_cost < not_visited[neighbor].cost)
                {
                    not_visited[neighbor].cost = new_cost;
                    not_visited[neighbor].parent = current;
                }
            }
        }
    }
    for (const auto& pair : not_visited)
    {
        shortest_path_tree[pair.first] = pair.second;
    }
}

void display_inital_shortest_path_tree(std::unordered_map<int, Couple>& visited)
{
    std::cout << "\nShortest Path tree from source vertex " << source_vertex << ":\n\n";
    for (const auto& pair : visited)
    {
        int vertex = pair.first;
        int cost = pair.second.cost;
        if (cost == INT_MAX)
        {
            std::cout << "Vertex " << vertex << " is unreachable from source vertex " << source_vertex << ".\n";
            continue;
        }
        std::vector<int> path;
        int current = vertex;
        while (current != -1)
        {
            path.push_back(current);
            current = visited[current].parent;
        }
        if (!path.empty())
        {
            std::cout << "Vertex: " << vertex << ", Cost: " << cost << ", Path: ";
            for (int i = path.size() - 1; i >= 0; --i)
            {
                std::cout << path[i];
                if (i != 0) std::cout << " -> ";
            }
            std::cout << '\n';
        }
    }
}

void insert_edge(int source, int destination, int* objectives, std::vector<Edge*>& edge_list)
{
    Edge* edge = new Edge();
    edge->source = source;
    edge->destination = destination;
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        edge->objectives[i] = objectives[i];
    }
    edge_list.push_back(edge);
    return;
}

void remove_edge(int source, int destination, std::vector<Edge*>& edge_list)
{
    for (unsigned int i = 0; i < edge_list.size(); ++i)
    {
        if (edge_list[i]->source == source && edge_list[i]->destination == destination)
        {
            delete edge_list[i];
            edge_list.erase(edge_list.begin() + i);
        }
    }
    return;
}

void display_edge_list(std::vector<Edge*>& edge_list)
{
    std::cout << '\n';
    for (Edge* edge : edge_list)
    {
        std::cout << edge->source << ", " << edge->destination << ", ";
        for (unsigned int j = 0; j < number_of_objectives; ++j)
        {
            std::cout << edge->objectives[j] << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    return;
}

void make_updated_shortest_path_tree(std::vector<Edge*>& edge_list, std::unordered_map<int, Couple>& shortes_path_tree)
{
    std::vector<Edge*> inserted_edges;

    return;
}

void display_updated_shortest_path_tree(std::unordered_map<int, Couple>& shortest_path_tree)
{
    return;
}

int main(int argc, char* argv[])
{
    std::string graph_file = read_graph_from_file("graph.txt");
    std::vector<Edge*> edge_list;
    create_edge_list(graph_file, edge_list);
    std::unordered_map<int, Couple> shortest_path_tree;
    make_inital_shortest_path_tree(edge_list, shortest_path_tree);
    display_inital_shortest_path_tree(shortest_path_tree);
    int* objectives = new int[1];
    *objectives = 8;
    std::cout << edge_list.size() << '\n';
    display_edge_list(edge_list);
    insert_edge(2, 4, objectives, edge_list);
    *objectives = 17;
    insert_edge(2, 1, objectives, edge_list);
    display_edge_list(edge_list);
    remove_edge(3, 4, edge_list);
    remove_edge(0, 1, edge_list);
    display_edge_list(edge_list);
    return 0;
}
