import pandas as pd 
import numpy as np 
import re
import os 
import glob
import pickle
import zipfile
import json
import ast
from fuzzywuzzy import fuzz
from datetime import datetime, timedelta
import torch
from transformers import AutoTokenizer, AutoModelForTokenClassification, AutoModelForCausalLM, BitsAndBytesConfig
import transformers
import torch
from difflib import SequenceMatcher

print('Starting...')

token = ''
from huggingface_hub import notebook_login
notebook_login()

print('After login')


# Load model directly
from transformers import AutoTokenizer, AutoModelForCausalLM


tokenizer = AutoTokenizer.from_pretrained("mistralai/Mistral-7B-Instruct-v0.2")
model = AutoModelForCausalLM.from_pretrained("mistralai/Mistral-7B-Instruct-v0.2", torch_dtype=torch.float16)

# Clear CUDA cache
torch.cuda.empty_cache()

# Check if a GPU is available and use it
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model.to(device)
 
# Print the device
print(f'Model is using device: {device}')

print('Done model download')

# Path to your folder containing the JSON files
folder_path = 'pathway_data/bioc'

# Lists to hold abstracts and their corresponding filenames
abstracts_list = []
filenames_list = []

# Paper sections to extract from
section_names = {"document title", "keywords", "Abstract"}
n = 0
# Loop through each file in the directory
for filename in os.listdir(folder_path):
    if filename.endswith('.json'):
        file_path = os.path.join(folder_path, filename)

        # Open and read the JSON file
        with open(file_path, 'r') as file:
            data = json.load(file)

        # Now, `data` is a dictionary containing the JSON data
        paper = data['documents'][0]
        text_data = paper['passages']

        # String to hold extracted biological_text values
        biological_text_list = ""

        # Loop through each dictionary in the list
        for entry in text_data:
            # Safely get the nested dictionary from the 'infons' key
            name_dict = entry.get('infons', {})
            
            # Check if the 'section_title_1' key's value is in the set of section_names
            if name_dict.get('section_title_1') in section_names:
                # Extract the value associated with 'text'
                biological_text = entry.get('text', '')

                # Check if the current string already ends with a full stop
                if not biological_text_list.endswith('. '):
                    biological_text_list += '. '
                
                # Append the biological_text value
                biological_text_list += biological_text

        # Ensure there is no leading full stop and space
        if biological_text_list.startswith('. '):
            biological_text_list = biological_text_list[2:]

        # Add a final full stop at the end if not present
        if not biological_text_list.endswith('.'):
            biological_text_list += '.'

        # Use re.sub() to replace multiple spaces with a single space
        biological_text_list = re.sub(r'\s+', ' ', biological_text_list)

        # Append the cleaned abstract to the list
        abstracts_list.append(biological_text_list)
        # Append the filename to the filenames list
        filenames_list.append(filename)

# Specify the directory to store the text files
output_dir = '../pathway_data/bioc_cancer_types'
os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist

# Function to get the set of base filenames (without extensions) from a directory
def get_base_filenames(folder_path, suffix='_result.txt'):
    base_filenames = set()
    for filename in os.listdir(folder_path):
        if filename.endswith(suffix):
            base_filename, _ = os.path.splitext(filename.replace(suffix, ''))
            base_filenames.add(base_filename)
    return base_filenames

# Get the set of base filenames from the output directory
processed_base_filenames = get_base_filenames(output_dir)
# Now we need to remove already processed files from abstracts_list and filenames_list
abstracts_list_filtered = []
filenames_list_filtered = []

for i in range(len(filenames_list)):
    base_filename, _ = os.path.splitext(filenames_list[i])
    if base_filename not in processed_base_filenames:
        abstracts_list_filtered.append(abstracts_list[i])
        filenames_list_filtered.append(filenames_list[i])

# Use the filtered lists for further processing
abstracts_list = abstracts_list_filtered
filenames_list = filenames_list_filtered
print(filenames_list)

print(f"Finished obtaining abstract text. {len(abstracts_list_filtered)} abstracts remaining to be processed.")

print("Finshed obtaining abstract text.")
print("Starting cancer classification query...")

# Empty list to contain the filename and cancer type pairs
file_cancer_list = []
results = []

# Regular expression to find the text between [/INST] and </s>
pattern = re.compile(r'\[/INST\](.*?)</s>', re.DOTALL)

# Maximum length for the prompt
max_prompt_length = 4000

# Function to split text into chunks ending at a full stop
def split_text_into_chunks(text, max_length):
    chunks = []
    while len(text) > max_length:
        split_pos = text.rfind('. ', 0, max_length)
        if split_pos == -1:
            split_pos = max_length
        chunk = text[:split_pos+1].strip()
        chunks.append(chunk)
        text = text[split_pos+1:].strip()
    if text:
        chunks.append(text)
    return chunks

# Loop through each abstract
for i in range(len(abstracts_list)):
    try:
        # Get the abstract
        abstract = abstracts_list[i]
        #print(abstract)
        # Split the abstract into chunks
        chunks = split_text_into_chunks(abstract, max_prompt_length)
        
        # Temporary list to store cancer types identified in each chunk
        temp_cancer_list = []

        for chunk in chunks:
            #print(chunk)
            chat = [
                {
                    "role": "system", 
                    "content": "You are highly intelligent Natural Language Processing (NLP) system, specialised in Biomedical Document Classification.",
                    "role": "user", 
                    "content": f"Based on the following abstract list me all the cancers that are mentioned and will be studied in the publications: {chunk}. You should output just list seperated by commas of the identified cancers with no other words or output 'There are no mention of cancers' if there are no mention of cancers."},
            ]
            
            prompt = tokenizer.apply_chat_template(chat, tokenize=False, add_generation_prompt=True)
            
            inputs = tokenizer.encode(prompt, add_special_tokens=False, return_tensors="pt")

            outputs = model.generate(input_ids=inputs.to(model.device), max_new_tokens=4096)
            print(tokenizer.decode(outputs[0]))
            
            print(f'Output tensors are on device: {outputs.device}')

            # Find all matches in the text
            matches = pattern.findall(tokenizer.decode(outputs[0]))

            # Loop through each match and check the condition
            for match in matches:
                extracted_text = match.strip().lower()

                if 'there are no mention of cancers' in extracted_text:
                    temp_cancer_list.append('There are no mention of cancers')
                else:
                    temp_cancer_list.append(match.strip())

        # Determine the final cancer type based on the temporary list
        if all('there are no mention of cancers' in entry.lower() for entry in temp_cancer_list):
            final_cancer_type = 'There are no mention of cancers'
        else:
            final_cancer_types = [entry for entry in temp_cancer_list if 'there are no mention of cancers' not in entry.lower()]
            final_cancer_type = ', '.join(final_cancer_types)
        
        print(final_cancer_type)

        # Split the filename to remove .json extension and add _result.txt
        base_filename, _ = os.path.splitext(filenames_list[i])
        text_filename = base_filename + '_result.txt'

        # Join the directory path with the filename
        output_path = os.path.join(output_dir, text_filename)

        # Write the result to a text file
        with open(output_path, 'w') as text_file:
            text_file.write(final_cancer_type)

        print(f"{i+1} out of {len(abstracts_list)} processed")
    
    except Exception as e:
        print(f"Error processing abstract {i+1}: {e}")
print("All files processed.")