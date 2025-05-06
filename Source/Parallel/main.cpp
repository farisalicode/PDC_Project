#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <climits>
#include <mpi.h>
#include <omp.h>



std::vector<std::unordered_map<int, Couple>> global_sosps;  
const unsigned int number_of_objectives = 8;
const int source_vertex = 0;

struct Edge
{
public:
    int source;
    int destination;
    int* objectives;
    Edge() {
        this->source = -1;
        this->destination = -1;
        this->objectives = new int[number_of_objectives];
    }
    Edge(int source, int destination) {
        this->source = source;
        this->destination = destination;
        this->objectives = new int[number_of_objectives];
    }
    Edge(const Edge& edge) {
        this->source = edge.source;
        this->destination = edge.destination;
        this->objectives = new int[number_of_objectives];
        for (unsigned int i = 0; i < number_of_objectives; ++i) {
            this->objectives[i] = edge.objectives[i];
        }
    }
    Edge& operator=(const Edge& edge) {
        this->source = edge.source;
        this->destination = edge.destination;
        for (unsigned int i = 0; i < number_of_objectives; ++i) {
            this->objectives[i] = edge.objectives[i];
        }
        return *this;
    }
    ~Edge() {
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
    Couple() {
        this->parent = -1;
        this->cost = -1;
    }
    Couple(int parent, int cost) {
        this->parent = parent;
        this->cost = cost;
    }
};

struct SerializedCouple
{
    int vertex;
    int parent;
    int cost;
};

std::string read_graph_from_file(const char* file_name)
{
    std::string graph_file;
    std::ifstream file(file_name);
    if (!file) {
        std::cerr << "Failed to open the file " << file_name << ".\n";
        return graph_file;
    }
    std::string line;
    while (std::getline(file, line)) {
        graph_file += line + '\n';
    }
    file.close();
    return graph_file;
}

void create_edge_list(const std::string& graph_file, std::vector<Edge*>& edge_list)
{
    std::string number = "";
    unsigned int comma_count = 0;
    unsigned int objective_count = 0;
    Edge* edge = new Edge();
    for (unsigned int i = 0; i < graph_file.length(); ++i) {
        if (graph_file[i] >= '0' && graph_file[i] <= '9') {
            number += graph_file[i];
        }
        else if (graph_file[i] == ',') {
            if (comma_count == 0) {
                edge->source = std::stoi(number);
                number = "";
                ++comma_count;
            }
            else if (comma_count == 1) {
                edge->destination = std::stoi(number);
                number = "";
                ++comma_count;
            }
            else if (comma_count > 1) {
                edge->objectives[objective_count] = std::stoi(number);
                number = "";
                ++objective_count;
            }
            continue;
        }
        else if (graph_file[i] == '\n') {
            edge_list.push_back(edge);
            comma_count = 0;
            objective_count = 0;
            edge = new Edge();
            continue;
        }
        else if (graph_file[i] == '~') {
            edge_list.push_back(edge);
            comma_count = 0;
            objective_count = 0;
            break;
        }
        else {
            std::cerr << "Invalid character found in graph file.\n";
            return;
        }
    }
    return;
}

void insert_edge(int source, int destination, int* objectives, std::vector<Edge*>& edge_list)
{
    Edge* edge = new Edge();
    edge->source = source;
    edge->destination = destination;
    for (unsigned int i = 0; i < number_of_objectives; ++i) {
        edge->objectives[i] = objectives[i];
    }
    edge_list.push_back(edge);
    return;
}

void remove_edge(int source, int destination, std::vector<Edge*>& edge_list)
{
    for (unsigned int i = 0; i < edge_list.size(); ++i) {
        if (edge_list[i]->source == source && edge_list[i]->destination == destination) {
            delete edge_list[i];
            edge_list.erase(edge_list.begin() + i);
        }
    }
    return;
}

void display_edge_list(std::vector<Edge*>& edge_list)
{
    std::cout << '\n';
    for (Edge* edge : edge_list) {
        std::cout << edge->source << ", " << edge->destination << ", ";
        for (unsigned int j = 0; j < number_of_objectives; ++j) {
            std::cout << edge->objectives[j] << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    return;
}

int find_minimum(std::unordered_map<int, Couple>& vertices)
{
    int minimum = INT_MAX;
    int vertex = -1;
    for (auto& pair : vertices) {
        if (pair.second.cost < minimum) {
            minimum = pair.second.cost;
            vertex = pair.first;
        }
    }
    return vertex;
}

void make_static_sosp(std::vector<Edge*>& edge_list, std::unordered_map<int, Couple>& sosp, unsigned int objective_number)
{
    std::unordered_map<int, Couple> not_visited;
    for (Edge* edge : edge_list) {
        not_visited[edge->source] = Couple(-1, INT_MAX);
        not_visited[edge->destination] = Couple(-1, INT_MAX);
    }
    not_visited[source_vertex].cost = 0;
    while (!not_visited.empty()) {
        int current = find_minimum(not_visited);
        if (current == -1) {
            break;
        }
        Couple current_data = not_visited[current];
        sosp[current] = current_data;
        not_visited.erase(current);
        if (current_data.cost == INT_MAX) {
            continue;
        }
        for (unsigned int i = 0; i < edge_list.size(); ++i) {
            if (edge_list[i]->source != current) {
                continue;
            }
            int neighbor = edge_list[i]->destination;
            if (not_visited.find(neighbor) != not_visited.end()) {
                int weight = edge_list[i]->objectives[objective_number];
                int new_cost = current_data.cost + weight;
                if (new_cost < not_visited[neighbor].cost) {
                    not_visited[neighbor].cost = new_cost;
                    not_visited[neighbor].parent = current;
                }
            }
        }
    }
    for (const auto pair : not_visited) {
        sosp[pair.first] = pair.second;
    }
    return;
}

void display_static_sosp(std::unordered_map<int, Couple>& sosp, unsigned int objective_number)
{
    std::cout << "\nShortest Path tree from source vertex " << source_vertex << ", objective number " << objective_number << ":\n\n";
    for (const auto pair : sosp) {
        int vertex = pair.first;
        int cost = pair.second.cost;
        if (cost == INT_MAX) {
            std::cout << "Vertex " << vertex << " is unreachable from source vertex " << source_vertex << ".\n";
            continue;
        }
        std::vector<int> path;
        int current = vertex;
        while (current != -1) {
            path.push_back(current);
            current = sosp[current].parent;
        }
        if (!path.empty()) {
            std::cout << "Vertex: " << vertex << ", Cost: " << cost << ", Path: ";
            for (int i = path.size() - 1; i >= 0; --i) {
                std::cout << path[i];
                if (i != 0) std::cout << " -> ";
            }
            std::cout << '\n';
        }
    }
    return;
}

void make_static_mosp(std::vector<Edge*>& edge_list, std::unordered_map<int, Couple>& mosp, int rank, int size)
{
    int local_obj_start = 0;
    int local_obj_count = 0;

    MPI_Datatype edge_type;
    int blocklengths[3] = { 1, 1, number_of_objectives };
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = { MPI_INT, MPI_INT, MPI_INT };

    Edge temp_edge;
    MPI_Aint base_address;
    MPI_Get_address(&temp_edge, &base_address);
    MPI_Get_address(&temp_edge.source, &displacements[0]);
    MPI_Get_address(&temp_edge.destination, &displacements[1]);
    MPI_Get_address(temp_edge.objectives, &displacements[2]);

    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);

    MPI_Type_create_struct(3, blocklengths, displacements, types, &edge_type);
    MPI_Type_commit(&edge_type);

    MPI_Datatype couple_type;
    MPI_Type_contiguous(3, MPI_INT, &couple_type);
    MPI_Type_commit(&couple_type);

    std::vector<std::unordered_map<int, Couple>> sosps;

    local_obj_count = number_of_objectives / size;
    if (rank < number_of_objectives % size) {
        local_obj_count++;
    }
    local_obj_start = (number_of_objectives / size) * rank;
    if (rank >= number_of_objectives % size) {
        local_obj_start += number_of_objectives % size;
    }
    else {
        local_obj_start += rank;
    }

    for (unsigned int i = 0; i < local_obj_count; ++i) {
        unsigned int objective_index = local_obj_start + i;
        if (objective_index < number_of_objectives) {
            sosps.emplace_back(std::unordered_map<int, Couple>());
            make_static_sosp(edge_list, sosps[i], objective_index);
        }
    }

    std::vector<SerializedCouple> serialized_couples;
    std::vector<int> couple_counts(size, 0);
    std::vector<int> couple_displacements(size, 0);
    int total_couples = 0;

    if (rank == 0) {
        for (int proc = 0; proc < size; proc++) {
            couple_counts[proc] = 0;
        }
    }

    int local_couple_count = 0;
    for (unsigned int i = 0; i < sosps.size(); ++i) {
        local_couple_count += sosps[i].size();
    }

    MPI_Gather(&local_couple_count, 1, MPI_INT, couple_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        couple_displacements[0] = 0;
        for (int i = 1; i < size; i++) {
            couple_displacements[i] = couple_displacements[i - 1] + couple_counts[i - 1];
        }
        total_couples = couple_displacements[size - 1] + couple_counts[size - 1];
        serialized_couples.resize(total_couples);
    }

    // Prepare local data for sending
    std::vector<SerializedCouple> local_couples;
    for (unsigned int i = 0; i < sosps.size(); ++i) {
        unsigned int objective_index = local_obj_start + i;
        for (const auto& pair : sosps[i]) {
            SerializedCouple sc;
            sc.vertex = pair.first;
            sc.parent = pair.second.parent;
            sc.cost = pair.second.cost;
            local_couples.push_back(sc);
        }
    }

    // Gather all serialized couples to root
    MPI_Gatherv(local_couples.data(), local_couples.size(), couple_type,
        serialized_couples.data(), couple_counts.data(), couple_displacements.data(),
        couple_type, 0, MPI_COMM_WORLD);

    // Root process combines the results
    if (rank == 0) {
        std::vector<std::unordered_map<int, Couple>> all_sosps(number_of_objectives);

        // Reconstruct all SOSPs from gathered data
        for (int i = 0; i < total_couples; ++i) {
            SerializedCouple& sc = serialized_couples[i];
            int objective_index = i % number_of_objectives;
            all_sosps[objective_index][sc.vertex] = Couple(sc.parent, sc.cost);
        }

        // Create the combined graph using the SOSP trees
        std::vector<Edge*> combined_graph;
        std::unordered_map<int, std::unordered_map<int, int>> edge_count;

        for (unsigned int i = 0; i < number_of_objectives; ++i) {
            for (const auto& pair : all_sosps[i]) {
                int child = pair.first;
                int parent = pair.second.parent;
                if (parent != -1) {
                    edge_count[parent][child]++;
                }
            }
        }

        for (const auto& source_pair : edge_count) {
            int source = source_pair.first;
            for (const auto& dest_pair : source_pair.second) {
                int destination = dest_pair.first;
                int x = dest_pair.second;
                Edge* new_edge = new Edge(source, destination);
                new_edge->objectives = new int[1];
                new_edge->objectives[0] = number_of_objectives - x + 1;
                combined_graph.push_back(new_edge);
            }
        }

        make_static_sosp(combined_graph, mosp, 0);
        display_static_sosp(mosp, 0);

        for (Edge* edge : combined_graph) {
            delete edge;
        }
    }

    MPI_Type_free(&edge_type);
    MPI_Type_free(&couple_type);
    return;
}

void make_dynamic_sosp(std::vector<Edge*>& edge_list, std::unordered_map<int, Couple>& sosp,
    unsigned int objective_number,
    const std::vector<Edge*>& changed_edges = std::vector<Edge*>(),
    const std::vector<bool>& is_insertion = std::vector<bool>())
{

    if (sosp.empty()) {
        make_static_sosp(edge_list, sosp, objective_number);
        return;
    }

    std::unordered_map<int, std::vector<std::pair<int, int>>> insert_groups;
    for (Edge* edge : edge_list) {
        int u = edge->source;
        int v = edge->destination;
        int weight = edge->objectives[objective_number];
        insert_groups[v].push_back(std::make_pair(u, weight));
    }

    std::vector<int> affected_vertices;

#pragma omp parallel
    {
        std::vector<int> thread_affected;

#pragma omp for
        for (int i = 0; i < static_cast<int>(insert_groups.size()); i++) {
            auto it = insert_groups.begin();
            std::advance(it, i);
            int v = it->first;
            bool updated = false;
            for (auto& edge_info : it->second) {
                int u = edge_info.first;
                int weight = edge_info.second;
                if (sosp.find(u) != sosp.end() && sosp[u].cost != INT_MAX) {
                    int new_cost = sosp[u].cost + weight;

#pragma omp critical
                    {
                        if (sosp.find(v) == sosp.end() || new_cost < sosp[v].cost) {
                            sosp[v].cost = new_cost;
                            sosp[v].parent = u;
                            updated = true;
                        }
                    }
                }
            }
            if (updated) {
                thread_affected.push_back(v);
            }
        }

#pragma omp critical
        {
            affected_vertices.insert(affected_vertices.end(), thread_affected.begin(), thread_affected.end());
        }
    }

    bool changes_made = true;
    while (changes_made && !affected_vertices.empty()) {
        changes_made = false;
        std::vector<int> neighbors;
        std::vector<bool> is_neighbor(1000, false);

        for (int affected : affected_vertices) {
            for (Edge* edge : edge_list) {
                if (edge->source == affected && !is_neighbor[edge->destination]) {
                    neighbors.push_back(edge->destination);
                    is_neighbor[edge->destination] = true;
                }
            }
        }

        affected_vertices.clear();

#pragma omp parallel
        {
            std::vector<int> thread_affected;
            bool thread_changes_made = false;

#pragma omp for
            for (int i = 0; i < static_cast<int>(neighbors.size()); i++) {
                int neighbor = neighbors[i];
                bool updated = false;
                for (Edge* edge : edge_list) {
                    if (edge->destination == neighbor && sosp.find(edge->source) != sosp.end()) {
                        int u = edge->source;
                        int weight = edge->objectives[objective_number];
                        if (sosp[u].cost != INT_MAX) {
                            int new_cost = sosp[u].cost + weight;

#pragma omp critical
                            {
                                if (sosp.find(neighbor) == sosp.end() || new_cost < sosp[neighbor].cost) {
                                    sosp[neighbor].cost = new_cost;
                                    sosp[neighbor].parent = u;
                                    updated = true;
                                }
                            }
                        }
                    }
                }
                if (updated) {
                    thread_affected.push_back(neighbor);
                    thread_changes_made = true;
                }
            }

#pragma omp critical
            {
                affected_vertices.insert(affected_vertices.end(), thread_affected.begin(), thread_affected.end());
                if (thread_changes_made) {
                    changes_made = true;
                }
            }
        }
    }
    return;
}

void display_dynamic_sosp(std::unordered_map<int, Couple>& sosp, unsigned int objective_number)
{
    std::cout << "\nShortest Path tree from source vertex " << source_vertex << ", objective number " << objective_number << ":\n\n";
    for (const auto pair : sosp) {
        int vertex = pair.first;
        int cost = pair.second.cost;
        if (cost == INT_MAX) {
            std::cout << "Vertex " << vertex << " is unreachable from source vertex " << source_vertex << ".\n";
            continue;
        }
        std::vector<int> path;
        int current = vertex;
        while (current != -1) {
            path.push_back(current);
            current = sosp[current].parent;
        }
        if (!path.empty()) {
            std::cout << "Vertex: " << vertex << ", Cost: " << cost << ", Path: ";
            for (int i = path.size() - 1; i >= 0; --i) {
                std::cout << path[i];
                if (i != 0) std::cout << " -> ";
            }
            std::cout << '\n';
        }
    }
    return;
}

void make_dynamic_mosp(std::vector<Edge*>& edge_list, std::unordered_map<int, Couple>& mosp,
    const std::vector<Edge*>& changed_edges, const std::vector<bool>& is_insertion,
    int rank, int size)
{

    int local_obj_start = 0;
    int local_obj_count = 0;

    MPI_Datatype couple_type;
    MPI_Type_contiguous(3, MPI_INT, &couple_type);
    MPI_Type_commit(&couple_type);

    std::vector<std::unordered_map<int, Couple>> sosps;

    local_obj_count = number_of_objectives / size;
    if (rank < number_of_objectives % size) {
        local_obj_count++;
    }
    local_obj_start = (number_of_objectives / size) * rank;
    if (rank >= number_of_objectives % size) {
        local_obj_start += number_of_objectives % size;
    }
    else {
        local_obj_start += rank;
    }

    for (unsigned int i = 0; i < local_obj_count; ++i) {
        unsigned int objective_index = local_obj_start + i;
        if (objective_index < number_of_objectives) {
            sosps.emplace_back(std::unordered_map<int, Couple>());
            make_dynamic_sosp(edge_list, sosps[i], objective_index, changed_edges, is_insertion);
        }
    }

    std::vector<SerializedCouple> serialized_couples;
    std::vector<int> couple_counts(size, 0);
    std::vector<int> couple_displacements(size, 0);
    int total_couples = 0;

    if (rank == 0) {
        for (int proc = 0; proc < size; proc++) {
            couple_counts[proc] = 0;
        }
    }

    int local_couple_count = 0;
    for (unsigned int i = 0; i < sosps.size(); ++i) {
        local_couple_count += sosps[i].size();
    }

    MPI_Gather(&local_couple_count, 1, MPI_INT, couple_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        couple_displacements[0] = 0;
        for (int i = 1; i < size; i++) {
            couple_displacements[i] = couple_displacements[i - 1] + couple_counts[i - 1];
        }
        total_couples = couple_displacements[size - 1] + couple_counts[size - 1];
        serialized_couples.resize(total_couples);
    }

    std::vector<SerializedCouple> local_couples;
    for (unsigned int i = 0; i < sosps.size(); ++i) {
        unsigned int objective_index = local_obj_start + i;
        for (const auto& pair : sosps[i]) {
            SerializedCouple sc;
            sc.vertex = pair.first;
            sc.parent = pair.second.parent;
            sc.cost = pair.second.cost;
            local_couples.push_back(sc);
        }
    }

    MPI_Gatherv(local_couples.data(), local_couples.size(), couple_type,
        serialized_couples.data(), couple_counts.data(), couple_displacements.data(),
        couple_type, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::vector<std::unordered_map<int, Couple>> all_sosps(number_of_objectives);

        for (int i = 0; i < total_couples; ++i) {
            SerializedCouple& sc = serialized_couples[i];
            int objective_index = i % number_of_objectives;
            all_sosps[objective_index][sc.vertex] = Couple(sc.parent, sc.cost);
        }

        std::vector<Edge*> combined_graph;
        std::unordered_map<int, std::unordered_map<int, int>> edge_count;

        for (unsigned int i = 0; i < number_of_objectives; ++i) {
            for (const auto& pair : all_sosps[i]) {
                int child = pair.first;
                int parent = pair.second.parent;
                if (parent != -1) {
                    edge_count[parent][child]++;
                }
            }
        }

        for (const auto& source_pair : edge_count) {
            int source = source_pair.first;
            for (const auto& dest_pair : source_pair.second) {
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

        for (Edge* edge : combined_graph) {
            delete edge;
        }
    }

    MPI_Type_free(&couple_type);
    return;
}

int main(int argc, char* argv[])
{
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string graph_file;
    std::vector<Edge*> edge_list;
    std::unordered_map<int, Couple> mosp;

    if (rank == 0) {
        graph_file = read_graph_from_file("graph.txt");
        create_edge_list(graph_file, edge_list);
    }

    int edge_count = 0;
    if (rank == 0) {
        edge_count = edge_list.size();
    }
    MPI_Bcast(&edge_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Datatype edge_type;
    int blocklengths[3] = { 1, 1, number_of_objectives };
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = { MPI_INT, MPI_INT, MPI_INT };

    Edge temp_edge;
    MPI_Aint base_address;
    MPI_Get_address(&temp_edge, &base_address);
    MPI_Get_address(&temp_edge.source, &displacements[0]);
    MPI_Get_address(&temp_edge.destination, &displacements[1]);
    MPI_Get_address(temp_edge.objectives, &displacements[2]);

    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);

    MPI_Type_create_struct(3, blocklengths, displacements, types, &edge_type);
    MPI_Type_commit(&edge_type);

    std::vector<int> edge_sources(edge_count);
    std::vector<int> edge_destinations(edge_count);
    std::vector<int> edge_objectives(edge_count * number_of_objectives);

    if (rank == 0) {
        for (int i = 0; i < edge_count; i++) {
            edge_sources[i] = edge_list[i]->source;
            edge_destinations[i] = edge_list[i]->destination;
            for (unsigned int j = 0; j < number_of_objectives; j++) {
                edge_objectives[i * number_of_objectives + j] = edge_list[i]->objectives[j];
            }
        }
    }

    MPI_Bcast(edge_sources.data(), edge_count, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(edge_destinations.data(), edge_count, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(edge_objectives.data(), edge_count * number_of_objectives, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        for (int i = 0; i < edge_count; i++) {
            Edge* edge = new Edge(edge_sources[i], edge_destinations[i]);
            for (unsigned int j = 0; j < number_of_objectives; j++) {
                edge->objectives[j] = edge_objectives[i * number_of_objectives + j];
            }
            edge_list.push_back(edge);
        }
    }

    make_static_mosp(edge_list, mosp, rank, size);

    std::vector<Edge*> changed_edges;
    std::vector<bool> is_insertion;

    if (rank == 0) {
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

        delete[] objectives;
    }

    int changed_edge_count = 0;
    if (rank == 0) {
        changed_edge_count = changed_edges.size();
    }
    MPI_Bcast(&changed_edge_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> changed_sources(changed_edge_count);
    std::vector<int> changed_destinations(changed_edge_count);
    std::vector<int> changed_objectives(changed_edge_count * number_of_objectives);
    std::vector<int> changed_is_insertion(changed_edge_count);

    if (rank == 0) {
        for (int i = 0; i < changed_edge_count; i++) {
            changed_sources[i] = changed_edges[i]->source;
            changed_destinations[i] = changed_edges[i]->destination;
            for (unsigned int j = 0; j < number_of_objectives; j++) {
                if (is_insertion[i]) {
                    changed_objectives[i * number_of_objectives + j] = changed_edges[i]->objectives[j];
                }
                else {
                    changed_objectives[i * number_of_objectives + j] = 0;
                }
            }
            changed_is_insertion[i] = is_insertion[i] ? 1 : 0;
        }
    }

    MPI_Bcast(changed_sources.data(), changed_edge_count, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(changed_destinations.data(), changed_edge_count, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(changed_objectives.data(), changed_edge_count * number_of_objectives, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(changed_is_insertion.data(), changed_edge_count, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        for (int i = 0; i < changed_edge_count; i++) {
            if (!changed_is_insertion[i]) {
                remove_edge(changed_sources[i], changed_destinations[i], edge_list);
            }
        }

        for (int i = 0; i < changed_edge_count; i++) {
            if (changed_is_insertion[i]) {
                int* objectives = new int[number_of_objectives];
                for (unsigned int j = 0; j < number_of_objectives; j++) {
                    objectives[j] = changed_objectives[i * number_of_objectives + j];
                }
                insert_edge(changed_sources[i], changed_destinations[i], objectives, edge_list);
                delete[] objectives;
            }
        }

        for (int i = 0; i < changed_edge_count; i++) {
            Edge* edge = new Edge(changed_sources[i], changed_destinations[i]);
            if (changed_is_insertion[i]) {
                for (unsigned int j = 0; j < number_of_objectives; j++) {
                    edge->objectives[j] = changed_objectives[i * number_of_objectives + j];
                }
            }
            changed_edges.push_back(edge);
            is_insertion.push_back(changed_is_insertion[i] == 1);
        }
    }

    make_dynamic_mosp(edge_list, mosp, changed_edges, is_insertion, rank, size);

    for (Edge* edge : changed_edges) {
        delete edge;
    }

    for (Edge* edge : edge_list) {
        delete edge;
    }

    MPI_Type_free(&edge_type);
    MPI_Finalize();
    return 0;
}
