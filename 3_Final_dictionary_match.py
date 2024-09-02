import pandas as pd
import spacy
from spacy.matcher import Matcher
from rapidfuzz import fuzz
import re
import glob
import os
import json
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed

# Load the SpaCy model once globally
nlp = spacy.load("en_core_web_sm")

# Define global variables for matchers and pattern lists
matcher_first = None
matcher_second = None
matcher_third = None
matcher_fourth = None
pattern_list_first = []
pattern_list_second = []
pattern_list_third = []
pattern_list_fourth = []


def initialize_matchers_and_patterns():
    global matcher_first, matcher_second, matcher_third, matcher_fourth
    global pattern_list_first, pattern_list_second, pattern_list_third, pattern_list_fourth

    df = pd.read_excel("./Final_metabolic_pathway_list_v4_with_ori_vstar.xlsx")
    
    # Function to check if a string contains brackets
    def has_brackets(name):
        return '(' in name and ')' in name

    # Create two separate dataframes based on the presence of brackets in 'Pathway_Name'
    with_brackets = df[df['Pathway_Name'].apply(has_brackets)]
    without_brackets = df[~df['Pathway_Name'].apply(has_brackets)]
    #print(with_brackets.head())
    #print(without_brackets.head())

    # Define exceptions
    exceptions = [
        'signalling pathway', 'metabolism pathway', 'biosynthesis pathway',
        'glycolysis pathway', 'degradation pathway', 'conversion pathway',
        'fixation pathway', 'elongation pathway', 'signaling pathway'
    ]

    # Function to create patterns for pathway names
    def create_patterns(pathway_name):
        # Remove content within parentheses and the parentheses themselves
        cleaned_name = re.sub(r'\s*\(.*?\)\s*', ' ', pathway_name).strip()
        tokens = cleaned_name.split()  # Split by spaces
        base_pattern = [{"LOWER": token.lower()} for token in tokens]
        patterns = [base_pattern]
        
        if tokens[-1].lower() == "pathway" and not any(cleaned_name.lower().endswith(exception) for exception in exceptions):
            patterns.extend([
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"LOWER": "and"}, {"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"IS_PUNCT": True}, {"IS_ALPHA": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"IS_PUNCT": True}, {"IS_ALPHA": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"IS_PUNCT": True}, {"IS_ALPHA": True}, {"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"LOWER": "pathway"}],
                [{"LOWER": token.lower()} for token in tokens[:-1]] + [{"IS_PUNCT": True}, {"IS_ALPHA": True}, {"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": "pathway"}]
            ])
        else:
            # Check if it ends with any of the exceptions
            for exception in exceptions:
                if cleaned_name.lower().endswith(exception):
                    exception_tokens = exception.split()
                    patterns.extend([
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"LOWER": "and"}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"IS_PUNCT": True}, {"IS_ALPHA": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"IS_PUNCT": True}, {"IS_ALPHA": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"IS_PUNCT": True}, {"IS_ALPHA": True}, {"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}],
                        [{"LOWER": token.lower()} for token in tokens[:-len(exception_tokens)]] + [{"IS_PUNCT": True}, {"IS_ALPHA": True}, {"IS_PUNCT": True}, {"LOWER": "and"}, {"IS_ALPHA": True}, {"IS_ALPHA": True}, {"LOWER": exception_tokens[0]}, {"LOWER": exception_tokens[1]}]
                    ])
                    break
        
        return patterns

    pattern_list_first = []
    pattern_list_second = []
    pattern_list_third = []
    pattern_list_fourth = []

    for index, row in without_brackets.iterrows():
        pathway_name = row['Pathway_Name']
        pattern_id = row['ID']
        if any(pathway_name.strip().endswith(exception) for exception in exceptions):
            patterns = create_patterns(pathway_name)
            for pattern in patterns:
                pattern_list_second.append({"ID": pattern_id, "pattern": pattern, "pathway_name": pathway_name})
        elif pathway_name.endswith('pathway') and not any(pathway_name.strip().endswith(exception) for exception in exceptions):
            patterns = create_patterns(pathway_name)
            for pattern in patterns:
                pattern_list_fourth.append({"ID": pattern_id, "pattern": pattern, "pathway_name": pathway_name})
        else:
            tokens = re.sub(r'\s*\(.*?\)\s*', ' ', pathway_name).strip().split()
            pattern = [{"LOWER": token} for token in tokens]
            pattern_list_third.append({"ID": pattern_id, "pattern": pattern, "pathway_name": pathway_name})

    for index, row in with_brackets.iterrows():
        pathway_name = row['Pathway_Name']
        pattern_id = row['ID']
        original_tokens = pathway_name.strip().split()
        original_pattern = [{"LOWER": token.lower()} for token in original_tokens]
        pattern_list_first.append({"ID": pattern_id, "pattern": original_pattern, "pathway_name": pathway_name})

    matcher_first = Matcher(nlp.vocab)
    matcher_second = Matcher(nlp.vocab)
    matcher_third = Matcher(nlp.vocab)
    matcher_fourth = Matcher(nlp.vocab)

    for item in pattern_list_first:
        matcher_first.add(item['ID'], [item['pattern']])

    for item in pattern_list_second:
        matcher_second.add(item['ID'], [item['pattern']])

    for item in pattern_list_third:
        matcher_third.add(item['ID'], [item['pattern']])

    for item in pattern_list_fourth:
        matcher_fourth.add(item['ID'], [item['pattern']])


# Function to perform fuzzy matching on a sentence
def fuzzy_match_sentence(sent, pathway_df, cleaned_positions, threshold=90):
    matched_pathways = []
    count = 0
    no_match_count = 0
    
    for token in sent:
        for index, row in pathway_df.iterrows():
            unique_id = row['ID']
            pathway = row['Pathway_Name']
            words = pathway.split()
            first_word = words[0]
            similarity = fuzz.ratio(token.text, first_word)
            if len(words) == 1:
                # For single-word pathways, require an exact match
                if similarity == 100:
                    match_end_idx = token.idx + len(token.text)
                    if match_end_idx < len(cleaned_positions):
                        extracted_text = sent.text[token.idx:match_end_idx]
                        print('found fuzzy')
                        matched_pathways.append((cleaned_positions[token.idx], cleaned_positions[match_end_idx], pathway, unique_id))
                        print(f"Extracted Text: {extracted_text}")
                        print(f"Pathway: {pathway}")
                        if pathway == None or pathway == 'None':
                            no_match_count += 1
            else:
                # For multi-word pathways, use the defined threshold
                if threshold <= similarity < 100:
                    match_text = token.text
                    match_end_idx = token.idx + len(token.text)
                    first_word_start_idx = token.idx
                    next_token_idx = token.i + 1
                    all_parts_matched = True
                    for part in pathway.split()[1:]:
                        part_matched = False
                        while next_token_idx < len(sent) and count < 4:
                            next_token = sent[next_token_idx]
                            match_text += ' ' + next_token.text
                            match_end_idx = next_token.idx + len(next_token.text)
                            next_token_idx += 1
                            count += 1
                            if fuzz.ratio(next_token.text, part) >= threshold:
                                part_matched = True
                                break
                        if not part_matched:
                            all_parts_matched = False
                            break
                    if all_parts_matched:
                        extracted_text = sent.text[token.idx:match_end_idx]
                        if '.' or ':' or ';' not in extracted_text[:-1]:
                            matched_pathways.append((cleaned_positions[token.idx], cleaned_positions[match_end_idx], pathway, unique_id))
                            print('found fuzzy')
                            print(f"Extracted Text: {extracted_text}")
                            print(f"Pathway: {pathway}")
                            if pathway == None or pathway == 'None':
                                no_match_count += 1
                        break  # Move to the next token once a match is found
    return matched_pathways, no_match_count


def process_sentence(sentence):
    # Create a list of positions for the original sentence
    original_positions = list(range(len(sentence)))

    # Use regex to find and remove parentheses and their content, along with the space after
    cleaned_sentence = re.sub(r'\(.*?\)\s*', '', sentence)

    # Find all matches of parentheses and their content
    matches = list(re.finditer(r'\(.*?\)\s*', sentence))

    cleaned_positions = []
    original_index = 0
    removed_positions = []

    # Record positions of removed content
    for match in matches:
        start, end = match.span()
        removed_positions.extend(range(start, end))

    # Adjust cleaned positions based on recorded removals
    for i, char in enumerate(sentence):
        if i not in removed_positions:
            cleaned_positions.append(i)

    return original_positions, cleaned_sentence, cleaned_positions

def find_matches(text, df, fuzzy_df):
    ori_text = text.lower()
    ori_text = ori_text + ' '
    #print(f"Original Text: {ori_text}")
    #print(ori_text)
    # Process the text
    doc = nlp(ori_text)
    no_fuzzy_match_count = 0

    # Initialize lists for matches
    exact_matches = []
    fuzzy_matches = []

    fuzzy_words = {'cycle', 'degradation', 'production', 'signalling', 'signaling', 'metabolism', 'pathway', 'synthesis', 'biosynthesis', 'cascade'}
    word_pattern = '|'.join(fuzzy_words)

    # Track the cumulative length of preceding sentences
    cumulative_length = 0

    # Split the text into sentences and process each sentence individually
    for sent in doc.sents:
        # Perform exact matching using matcher_first first
        matches_first = matcher_first(sent)
        for match_id, start, end in matches_first:
            span = sent[start:end]
            pattern_id = nlp.vocab.strings[match_id]
            pattern = next(item['pattern'] for item in pattern_list_first if item['ID'] == pattern_id)
            print(ori_text[span.start_char:span.end_char])
            exact_matches.append((span.start_char, span.end_char, 'None', nlp.vocab.strings[match_id]))
            print('Matcher first')

        # Update the cumulative length for the next sentence
        cumulative_length += len(sent.text)   # +1 accounts for the space or punctuation separating sentences

    # Re-process the text for cleaned positions
    original_positions, cleaned_text, cleaned_positions = process_sentence(ori_text)
    doc = nlp(cleaned_text)
    cumulative_length = 0  # Reset cumulative length for cleaned text
    for sent in doc.sents:
        # Perform exact matching using matcher_second first
        matches_second = matcher_second(sent)
        for match_id, start, end in matches_second:
            span = sent[start:end]
            pattern_id = nlp.vocab.strings[match_id]
            pattern = next(item['pattern'] for item in pattern_list_second if item['ID'] == pattern_id)
            print(ori_text[cleaned_positions[span.start_char]:cleaned_positions[span.end_char]])
            exact_matches.append((cleaned_positions[span.start_char], cleaned_positions[span.end_char], 'None', nlp.vocab.strings[match_id]))
            print('Matcher second')

        # Perform exact matching using matcher_third
        matches_third = matcher_third(sent)
        for match_id, start, end in matches_third:
            span = sent[start:end]
            pattern_id = nlp.vocab.strings[match_id]
            pattern = next(item['pattern'] for item in pattern_list_third if item['ID'] == pattern_id)
            print(ori_text[cleaned_positions[span.start_char]:cleaned_positions[span.end_char]])
            exact_matches.append((cleaned_positions[span.start_char], cleaned_positions[span.end_char], 'None', nlp.vocab.strings[match_id]))
            print('Matcher third')

        matches_fourth = matcher_fourth(sent)
        for match_id, start, end in matches_fourth:
            span = sent[start:end]
            pattern_id = nlp.vocab.strings[match_id]
            pattern = next(item['pattern'] for item in pattern_list_fourth if item['ID'] == pattern_id)
            exact_matches.append((cleaned_positions[span.start_char], cleaned_positions[span.end_char], 'None', nlp.vocab.strings[match_id]))
            print(ori_text[cleaned_positions[span.start_char]:cleaned_positions[span.end_char]])
            #except IndexError as e:
                #exact_matches.append((cleaned_positions[span.start_char], cleaned_positions[span.end_char-1], 'None', nlp.vocab.strings[match_id]))
                #print(ori_text[cleaned_positions[span.start_char]:cleaned_positions[span.end_char-1]])
            #exact_matches.append((cleaned_positions[span.start_char], cleaned_positions[span.end_char], 'None', nlp.vocab.strings[match_id]))
            print('Matcher fourth') 

        # Update the cumulative length for the next sentence
        #cumulative_length += len(sent.text)   # +1 accounts for the space or punctuation separating sentences
        #print(cumulative_length)
        if any(word in sent.text for word in fuzzy_words):
            # Get the list of IDs in the exact matches list
            exact_match_ids = [match[3] for match in exact_matches]

            # Create a copy of the ID column from the fuzzy match dataframe and remove all *
            fuzzy_df_ids_cleaned = fuzzy_df['ID'].str.replace('*', '', regex=False)

            # Create a copy of the exact_match_ids list and remove all *
            exact_match_ids_cleaned = [id_.replace('*', '') for id_ in exact_match_ids]

            # Remove these cleaned IDs and their associated pathway names from the pathway_dict
            filtered_pathway_df = fuzzy_df[~fuzzy_df_ids_cleaned.isin(exact_match_ids_cleaned)]
            # Remove these IDs and their associated pathway names from the pathway_dict
            #filtered_pathway_df = fuzzy_df[~fuzzy_df['ID'].isin(exact_match_ids)]

            # Filter rows in `filtered_pathway_df` to only include those with the specified words
            final_filtered_df = filtered_pathway_df[filtered_pathway_df['Pathway_Name'].str.contains(word_pattern, case=False, regex=True)]

            # Perform fuzzy matching
            if not final_filtered_df.empty:
                matched, no_matches = fuzzy_match_sentence(sent, final_filtered_df, cleaned_positions)
                fuzzy_matches.extend(matched)
                no_fuzzy_match_count += no_matches

    # Combine exact matches and fuzzy matches without removing duplicates
    combined_matches = exact_matches + fuzzy_matches
    return combined_matches, no_fuzzy_match_count

def process_file(file_path, formatted_date, pathway_df, fuzzy_df):
    initialize_matchers_and_patterns()
    try:
        #initialize_matchers_and_patterns()
        total = 1
        count = 0
        current_time = datetime.now() 
        formatted_date = current_time.strftime("%Y-%m-%dT%H:%M:%SZ")
        total_no_fuzzy_matches = 0
        with open(file_path, 'r', encoding='utf-8') as f:
            all_data = json.load(f)
        for j in range(len(all_data['documents'][0]['passages'])):
            annotation = []
            matches, no_fuzzy_matches = find_matches(all_data['documents'][0]['passages'][j]['text'], pathway_df, fuzzy_df)
            passage_text = all_data['documents'][0]['passages'][j]['text']
            passage = all_data['documents'][0]['passages'][j]
            total_no_fuzzy_matches += no_fuzzy_matches
            #print(matches)
            for start, end, word, identifier in matches:
                annotation.append([start, end, word, 'pathway', identifier])
            annotation_final = []
            for k in range(len(annotation)):    
                        #print('found annotation')
                        offset = annotation[k][0] + all_data['documents'][0]['passages'][j]['offset']
                        length = annotation[k][1] - annotation[k][0]
                        extracted_text = passage_text[offset - passage['offset']:offset - passage['offset'] + length]
                        id_count = len(annotation[k][4].split(', '))
                        if id_count < 11:
                            identifier = annotation[k][4]
                        else:
                            identifier = extracted_text
                        annotation_final.append({
                            'id' : str(total),
                            'database_id' : annotation[k][4],
                            'infons': {
                                'type': 'pathway',
                                'identifier': identifier,
                                'annotator': 'pathway@codiet.eu',
                                'updated_at': f'{formatted_date}'},
                            'text': annotation[k][2],
                            'locations': [{'offset': annotation[k][0] + all_data['documents'][0]['passages'][j]['offset'], 'length': annotation[k][1] - annotation[k][0]}]
                        })
                        total += 1
            all_data['documents'][0]['passages'][j]['annotations'] = annotation_final

        #output_path = f'./dict_output_v9/{os.path.basename(file_path)}'
        output_path = f'./dict_output_J_v1/{os.path.basename(file_path)}'
        #output_path = f'./test_data_annotated_redo/batch_1_annotated_final_code/{os.path.basename(file_path)}'
        print(f"Saving file to {output_path}")

        # Save the updated JSON file
        with open(output_path, 'w', encoding="utf-8") as fp:
            json.dump(all_data, fp, indent=2, ensure_ascii=False)

        # Update 'None' annotations
        updated = False
        for passage in all_data['documents'][0]['passages']:
            passage_text = passage['text']
            if 'annotations' in passage and passage['annotations']:
                for annotation in passage['annotations']:
                    if annotation['text'] is None or annotation['text'].lower() == 'none':
                        database_id = annotation['database_id']
                        possible_pathways = pathway_df[pathway_df['ID'] == database_id]['Pathway_Name']
                        if not possible_pathways.empty:
                            doc = nlp(passage['text'])
                            offset = annotation['locations'][0]['offset']
                            length = annotation['locations'][0]['length']
                            extracted_text_ori = passage_text[offset - passage['offset']:offset - passage['offset'] + length]
                            print(f"Extracted Pathway: {extracted_text_ori}")

                            # Find the best match among the possible pathways
                            best_match = None
                            highest_similarity = 0
                            for pathway in possible_pathways:
                                pathway = pathway.lower()
                                extracted_text = extracted_text_ori.lower()
                                words = pathway.split()
                                first_word = words[0]
                                similarity = fuzz.ratio(extracted_text.split()[0], first_word)
                                if len(words) == 1:
                                    # For single-word pathways, require an exact match
                                    if similarity >= 90:
                                        match_end_idx = len(extracted_text.split()[0])
                                        #best_match = pathway
                                        best_match = extracted_text_ori
                                        highest_similarity = similarity
                                        break
                                else:
                                    # For multi-word pathways, use the defined threshold
                                    if similarity >= 80:
                                        match_text = extracted_text.split()[0]
                                        match_end_idx = len(match_text)
                                        next_token_idx = 1
                                        all_parts_matched = True
                                        for part in pathway.split()[1:]:
                                            part_matched = False
                                            while next_token_idx < len(extracted_text.split()):
                                                next_token = extracted_text.split()[next_token_idx]
                                                match_text += ' ' + next_token
                                                match_end_idx += len(next_token) + 1  # +1 for the space
                                                next_token_idx += 1
                                                if fuzz.ratio(next_token, part) >= 80:
                                                    part_matched = True
                                                    break
                                            if not part_matched:
                                                all_parts_matched = False
                                                break
                                        if all_parts_matched:
                                            similarity = fuzz.ratio(extracted_text, pathway)
                                            if similarity > highest_similarity:
                                                highest_similarity = similarity
                                                #best_match = pathway
                                                best_match = extracted_text_ori

                            if best_match:
                                print(f"Best Match: {best_match} (Similarity: {highest_similarity})")
                                annotation['text'] = best_match
                                updated = True

        if updated:
            # Save the updated JSON file
            with open(output_path, 'w', encoding='utf-8') as fp:
                json.dump(all_data, fp, indent=2, ensure_ascii=False)
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return 0  # Ensure an integer is returned even if an error occurs

def main(num_workers=None, batch_size=1000):
    total_no_fuzzy_matches = 0

    start_time = datetime.now()
    print('Starting')

    #input_list = glob.glob('./trial_new_dictionary_match/*.json') 
    
    ##input_list = glob.glob('./pathway_data_complete/bioc/*.json')
    
    #input_list = glob.glob('./test_data_redo/batch_1/*.json')
    # Only keep the first 10 files
    #input_list = input_list[:10]
    # Calculate the midpoint of the list
    #midpoint = len(input_list) // 2  # Floor division to get an integer

    # Select the first half of the files
    #input_list = input_list[:midpoint]
    # Select the second half of the files
    #input_list = input_list[midpoint:]

    input_folder = './pathway_data_complete/bioc/'
    output_folder = './dict_output_J_v1/'

    # Get a list of all JSON files in the input folder
    input_list_folder = glob.glob(os.path.join(input_folder, '*.json'))

    # Get a list of all JSON files in the output folder
    output_list_folder = glob.glob(os.path.join(output_folder, '*.json'))

    # Get the filenames without the directory path
    input_filenames = [os.path.basename(file) for file in input_list_folder]
    output_filenames = [os.path.basename(file) for file in output_list_folder]

    # Find the files that are in the input folder but not in the output folder
    files_to_process = [file for file in input_filenames if file not in output_filenames]

    # Add the directory path back to the filenames
    input_list = [os.path.join(input_folder, file) for file in files_to_process]

    # Now, files_to_process contains only the files that need to be processed
    print(f"Files to process: {len(input_list)}")

    formatted_date = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    df = pd.read_excel("./Final_metabolic_pathway_list_v4_with_ori_vstar.xlsx")

    # Function to check if a string contains brackets
    def has_no_brackets(name):
        return '(' not in name and ')' not in name
    # Function to check if the number of words is between 2 and 4
    def word_count_in_range(name, min_words=2, max_words=2):
        word_count = len(name.split())
        return min_words <= word_count <= max_words

    # Filter the DataFrame
    fuzzy_df = df[df['Pathway_Name'].apply(has_no_brackets) & df['Pathway_Name'].apply(word_count_in_range)]
    #fuzzy_df = df

    # Divide input list into smaller batches
    batches = [input_list[i:i + batch_size] for i in range(0, len(input_list), batch_size)]

    for batch_idx, batch in enumerate(batches):
        print(f'Processing batch {batch_idx + 1} of {len(batches)}')
        
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_file, file_path, formatted_date, df, fuzzy_df): i for i, file_path in enumerate(batch)}
            for future in as_completed(futures):
                i = futures[future]
                try:
                    #result = future.result()
                    print(f'Processing batch {batch_idx + 1} of {len(batches)}')
                    print(f'Completed processing index: {i+1} of {len(batch)}')
                    #total_no_fuzzy_matches += result
                except Exception as e:
                    print(f"Error in future result for index {i+1}: {e}")

    end_time = datetime.now()
    elapsed_time = end_time - start_time
    print(f"Elapsed Time: {elapsed_time}")
    print(f"Total No Fuzzy Matches: {total_no_fuzzy_matches}")
    print('Finished')

if __name__ == "__main__":
    num_workers = 21  # Set the desired number of CPU cores
    main(num_workers=num_workers)

