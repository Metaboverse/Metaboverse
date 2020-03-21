import sys
import json
import re
import os

def process_data(data_path):
    with open(data_path) as json_file:
        net_data = json.load(json_file)
    nodes = net_data['nodes']
    links = net_data['links']
    pathway_dict = net_data['pathway_dictionary']
    super_pathways = net_data['super_pathways']
    reactions_dict = net_data['reaction_dictionary']

    reaction_nodes = []
    nodes_dict = {}

    # expression_values = []

    for node in nodes:
        # if node['values'] not in expression_values:
        #     expression_values.append(node['values'])
        if 'notes' in node.keys():
            del node['notes']
        nodes_dict[node['id']] = node 
        if node['type'] == 'reaction':
            reaction_nodes.append(node)
    
    ### only keep reactant & product links ###       
    links_new = []
    for link in links:
        if link['type'] in ['reactant', 'product']:
            links_new.append(link)

    ### assign links to reaction nodes ###
    reaction_nodes_dict = {}
    for node in reaction_nodes:
        node['links'] = {'reactant':[], 'product':[]}
        node['pathways'] = []
        reaction_nodes_dict[node['id']] = node
    for link in links_new:
        if link['source'] in reaction_nodes_dict.keys():
            reaction_nodes_dict[link['source']]['links']['product'].append(link)
        elif link['target'] in reaction_nodes_dict.keys():
            reaction_nodes_dict[link['target']]['links']['reactant'].append(link)   
    ### assign pathways to reaction nodes ###
    for p_key in pathway_dict:
        for node_id in pathway_dict[p_key]['reactions']:
            if node_id in reaction_nodes_dict.keys():
                reaction_nodes_dict[node_id]['pathways'].append(p_key)
            else:
                pathway_dict[p_key]['reactions'].remove(node_id)

    ### other nodes: only reactant & product nodes ###
    other_nodes_dict = {}
    for node_id in reaction_nodes_dict:
        react_node = reaction_nodes_dict[node_id]
        for link in react_node['links']['reactant']:
            other_node = nodes_dict[link['source']]
            if type(other_node['values']) == list:
                other_node['expression'] = other_node['values'][0]
            else:
                other_node['expression'] = "None"
            other_nodes_dict[link['source']] = other_node
        for link in react_node['links']['product']:
            other_node = nodes_dict[link['target']]
            if type(other_node['values']) == list:
                other_node['expression'] = other_node['values'][0]
            else:
                other_node['expression'] = "None"
            other_nodes_dict[link['target']] = other_node  

      
    
    net_data_new = {"react_nodes":reaction_nodes_dict, "other_nodes":other_nodes_dict, "pathways":pathway_dict}    
    with open(os.getcwd()+'/data/processed_data.json','w') as f:
        json.dump(net_data_new, f)
    

def main():
    path = sys.argv[1]
    process_data(path)
    print("finish data processing")

if __name__ == '__main__':
    main()