#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>



const unsigned int number_of_objectives = 2;
const int source_vertex = 0;

struct Edge
{
public:
    int source;
    int destination;
    int* objectives;
    Edge ()
    {
        this->source = -1;
        this->destination = -1;
        this->objectives = new int[number_of_objectives];
    }
    Edge (int source, int destination)
    {
        this->source = source;
        this->destination = destination;
        this->objectives = new int[number_of_objectives];
    }
    Edge (const Edge& edge)
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
    ~Edge ()
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
    Couple ()
    {
        this->parent = -1;
        this->cost = -1;
    }
    Couple (int parent, int cost)
    {
        this->parent = parent;
        this->cost = cost;
    }
};

std::string read_graph_from_file (const char* file_name)
{
    std::string graph_file;
    std::ifstream file (file_name);
    if (!file)
    {
        std::cerr << "Failed to open the file " << file_name << ".\n";
        return graph_file;
    }
    std::string line;
    while (std::getline (file, line))
    {
        graph_file += line + '\n';
    }
    file.close ();
    return graph_file;
}

void create_edge_list (const std::string& graph_file, std::vector<Edge*>& edge_list)
{
    std::string number = "";
    unsigned int comma_count = 0;
    unsigned int objective_count = 0;
    Edge* edge = new Edge ();
    for (unsigned int i = 0; i < graph_file.length (); ++i)
    {
        if (graph_file[i] >= '0' && graph_file[i] <= '9')
        {
            number += graph_file[i];
        }
        else if (graph_file[i] == ',')
        {
            if (comma_count == 0)
            {
                edge->source = std::stoi (number);
                number = "";
                ++comma_count;
            }
            else if (comma_count == 1)
            {
                edge->destination = std::stoi (number);
                number = "";
                ++comma_count;
            }
            else if (comma_count > 1)
            {
                edge->objectives[objective_count] = std::stoi (number);
                number = "";
                ++objective_count;
            }
            continue;
        }
        else if (graph_file[i] == '\n')
        {
            edge_list.push_back (edge);
            comma_count = 0;
            objective_count = 0;
            edge = new Edge ();
            continue;
        }
        else if (graph_file[i] == '~')
        {
            edge_list.push_back (edge);
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

void insert_edge (int source, int destination, int* objectives, std::vector<Edge*>& edge_list)
{
    Edge* edge = new Edge ();
    edge->source = source;
    edge->destination = destination;
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        edge->objectives[i] = objectives[i];
    }
    edge_list.push_back (edge);
    return;
}

void remove_edge(int source, int destination, std::vector<Edge*>& edge_list)
{
    for (unsigned int i = 0; i < edge_list.size (); ++i)
    {
        if (edge_list[i]->source == source && edge_list[i]->destination == destination)
        {
            delete edge_list[i];
            edge_list.erase (edge_list.begin() + i);
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

int find_minimum (std::unordered_map<int, Couple>& vertices)
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

void make_static_sosp (std::vector<Edge*>& edge_list, std::unordered_map<int, Couple>& sosp, unsigned int objective_number)
{
    std::unordered_map<int, Couple> not_visited;
    for (Edge* edge : edge_list)
    {
        not_visited[edge->source] = Couple (-1, INT_MAX);
        not_visited[edge->destination] = Couple (-1, INT_MAX);
    }
    not_visited[source_vertex].cost = 0;
    while (!not_visited.empty ())
    {
        int current = find_minimum (not_visited);
        if (current == -1)
        {
            break;
        }
        Couple current_data = not_visited[current];
        sosp[current] = current_data;
        not_visited.erase (current);
        if (current_data.cost == INT_MAX)
        {
            continue;
        }
        for (unsigned int i = 0; i < edge_list.size (); ++i)
        {
            if (edge_list[i]->source != current)
            {
                continue;
            }
            int neighbor = edge_list[i]->destination;
            if (not_visited.find (neighbor) != not_visited.end ())
            {
                int weight = edge_list[i]->objectives[objective_number];
                int new_cost = current_data.cost + weight;
                if (new_cost < not_visited[neighbor].cost)
                {
                    not_visited[neighbor].cost = new_cost;
                    not_visited[neighbor].parent = current;
                }
            }
        }
    }
    for (const auto pair : not_visited)
    {
        sosp[pair.first] = pair.second;
    }
    return;
}

void display_static_sosp (std::unordered_map<int, Couple>& sosp, unsigned int objective_number)
{
    std::cout << "\nShortest Path tree from source vertex " << source_vertex << ", objective number " << objective_number << ":\n\n";
    for (const auto pair : sosp)
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
            path.push_back (current);
            current = sosp[current].parent;
        }
        if (!path.empty ())
        {
            std::cout << "Vertex: " << vertex << ", Cost: " << cost << ", Path: ";
            for (int i = path.size () - 1; i >= 0; --i)
            {
                std::cout << path[i];
                if (i != 0) std::cout << " -> ";
            }
            std::cout << '\n';
        }
    }
    return;
}

void make_static_mosp (std::vector<Edge*>& edge_list, std::unordered_map<int, Couple> &mosp)
{
    std::vector<std::unordered_map<int, Couple>> sosps;
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        sosps.emplace_back (std::unordered_map<int, Couple>());
        make_static_sosp (edge_list, sosps[i], i);
    }
    std::vector<Edge*> combined_graph;
    std::unordered_map<int, std::unordered_map<int, int>> edge_count;
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        for (const auto& pair : sosps[i])
        {
            int child = pair.first;
            int parent = pair.second.parent;
            if (parent != -1)
            {
                edge_count[parent][child]++;
            }
        }
    }
    for (const auto& source_pair : edge_count)
    {
        int source = source_pair.first;
        for (const auto& dest_pair : source_pair.second)
        {
            int destination = dest_pair.first;
            int x = dest_pair.second;
            Edge* new_edge = new Edge(source, destination);
            new_edge->objectives = new int[1];
            new_edge->objectives[0] = number_of_objectives - x + 1;
            combined_graph.push_back (new_edge);
        }
    }
    make_static_sosp (combined_graph, mosp, 0);
    display_static_sosp (mosp, 0);
    for (Edge* edge : combined_graph)
    {
        delete edge;
    }
    return;
}

void make_dynamic_sosp(std::vector<Edge*>& edge_list, std::unordered_map<int, Couple>& sosp,
    unsigned int objective_number,
    const std::vector<Edge*>& changed_edges = std::vector<Edge*>(),
    const std::vector<bool>& is_insertion = std::vector<bool>())
{
    if (sosp.empty())
    {
        make_static_sosp(edge_list, sosp, objective_number);
        return;
    }
    // This part can be easily parallelized with OpenMP later
    //{
        // Begin parallel section - each group of vertices will be processed by a separate thread
        std::unordered_map<int, std::vector<std::pair<int, int>>> insert_groups;
        for (Edge* edge : edge_list)
        {
            int u = edge->source;
            int v = edge->destination;
            int weight = edge->objectives[objective_number];
            insert_groups[v].push_back(std::make_pair(u, weight));
        }

        std::vector<int> affected_vertices;
        for (auto& group : insert_groups)
        {
            int v = group.first;
            bool updated = false;
            for (auto& edge_info : group.second)
            {
                int u = edge_info.first;
                int weight = edge_info.second;
                if (sosp.find(u) != sosp.end() && sosp[u].cost != INT_MAX)
                {
                    int new_cost = sosp[u].cost + weight;
                    if (sosp.find(v) == sosp.end() || new_cost < sosp[v].cost)
                    {
                        sosp[v].cost = new_cost;
                        sosp[v].parent = u;
                        updated = true;
                    }
                }
            }
            if (updated)
            {
                affected_vertices.push_back(v);
            }
        //}
        // End parallel section
    }

    // Propagate updates through the graph
    bool changes_made = true;
    while (changes_made && !affected_vertices.empty())
    {
        changes_made = false;
        std::vector<int> neighbors;
        std::vector<bool> is_neighbor(1000, false);

        for (int affected : affected_vertices)
        {
            for (Edge* edge : edge_list)
            {
                if (edge->source == affected && !is_neighbor[edge->destination])
                {
                    neighbors.push_back(edge->destination);
                    is_neighbor[edge->destination] = true;
                }
            }
        }

        affected_vertices.clear();
        for (int neighbor : neighbors)
        {
            bool updated = false;
            for (Edge* edge : edge_list)
            {
                if (edge->destination == neighbor && sosp.find(edge->source) != sosp.end())
                {
                    int u = edge->source;
                    int weight = edge->objectives[objective_number];
                    if (sosp[u].cost != INT_MAX)
                    {
                        int new_cost = sosp[u].cost + weight;
                        if (sosp.find(neighbor) == sosp.end() || new_cost < sosp[neighbor].cost)
                        {
                            sosp[neighbor].cost = new_cost;
                            sosp[neighbor].parent = u;
                            updated = true;
                            changes_made = true;
                        }
                    }
                }
            }
            if (updated)
            {
                affected_vertices.push_back(neighbor);
            }
        }
    }
    return;
}

void display_dynamic_sosp(std::unordered_map<int, Couple>& sosp, unsigned int objective_number)
{
    std::cout << "\nShortest Path tree from source vertex " << source_vertex << ", objective number " << objective_number << ":\n\n";
    for (const auto pair : sosp)
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
            current = sosp[current].parent;
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
    return;
}

void make_dynamic_mosp(std::vector<Edge*>& edge_list, std::unordered_map<int, Couple>& mosp,
    const std::vector<Edge*>& changed_edges, const std::vector<bool>& is_insertion)
{
    std::vector<std::unordered_map<int, Couple>> sosps(number_of_objectives);

    // Update each SOSP tree for each objective
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        make_dynamic_sosp(edge_list, sosps[i], i, changed_edges, is_insertion);
        // display_dynamic_sosp(sosps[i], i);
    }

    // Create a combined graph using the updated SOSP trees
    std::vector<Edge*> combined_graph;
    std::unordered_map<int, std::unordered_map<int, int>> edge_count;

    // Count how many times each edge appears in SOSP trees
    for (unsigned int i = 0; i < number_of_objectives; ++i)
    {
        for (const auto& pair : sosps[i])
        {
            int child = pair.first;
            int parent = pair.second.parent;
            if (parent != -1)
            {
                edge_count[parent][child]++;
            }
        }
    }

    // Create the combined graph with weights based on edge appearances in SOSP trees
    for (const auto& source_pair : edge_count)
    {
        int source = source_pair.first;
        for (const auto& dest_pair : source_pair.second)
        {
            int destination = dest_pair.first;
            int x = dest_pair.second;
            Edge* new_edge = new Edge(source, destination);
            new_edge->objectives = new int[1];
            new_edge->objectives[0] = number_of_objectives - x + 1;
            combined_graph.push_back(new_edge);
        }
    }
    mosp.clear();
    make_static_sosp(combined_graph, mosp, 0);
    display_dynamic_sosp(mosp, 0);
    for (Edge* edge : combined_graph)
    {
        delete edge;
    }
    return;
}

int main(int argc, char* argv[])
{
    std::string graph_file = read_graph_from_file("graph.txt");
    std::vector<Edge*> edge_list;
    create_edge_list(graph_file, edge_list);
    std::unordered_map<int, Couple> mosp;
    make_static_mosp(edge_list, mosp);
    std::vector<Edge*> changed_edges;
    std::vector<bool> is_insertion;
    int* objectives = new int[2];
    objectives[0] = 17;
    objectives[1] = 45;
    Edge* new_edge1 = new Edge(2, 4);
    new_edge1->objectives[0] = objectives[0];
    new_edge1->objectives[1] = objectives[1];
    changed_edges.push_back(new_edge1);
    is_insertion.push_back(true);
    insert_edge(2, 4, objectives, edge_list);
    Edge* new_edge2 = new Edge(2, 1);
    new_edge2->objectives[0] = objectives[0];
    new_edge2->objectives[1] = objectives[1];
    changed_edges.push_back(new_edge2);
    is_insertion.push_back(true);
    insert_edge(2, 1, objectives, edge_list);
    Edge* removed_edge1 = new Edge(3, 4);
    changed_edges.push_back(removed_edge1);
    is_insertion.push_back(false);
    remove_edge(3, 4, edge_list);
    Edge* removed_edge2 = new Edge(0, 1);
    changed_edges.push_back(removed_edge2);
    is_insertion.push_back(false);
    remove_edge(0, 1, edge_list);
    make_dynamic_mosp(edge_list, mosp, changed_edges, is_insertion);
    delete[] objectives;
    for (Edge* edge : changed_edges)
    {
        delete edge;
    }
    for (Edge* edge : edge_list)
    {
        delete edge;
    }
    return 0;
}
