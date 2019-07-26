#!/usr/bin/env python
import sys

#this is the annotation class. Annotations have a feature id, a key, and a value (in that order).
class Annotation:

    def __init__(self, feature_id='', key='', value=''):
        self.feature_id = feature_id
        self.key = key
        self.value = value

    #annotations are equal if their feature id's, keys, and values are equal
    def __eq__(self, other):
        if self.feature_id == other.feature_id and\
           self.key == other.key and\
           self.value == other.value:
            return True
        return False

    #they aren't equal...if they are not equal
    def __ne__(self, other):
        return not self == other

    #we order an annotation by first comparing feature id's, then by comparing keys if the feature id's are equal, then by values if the keys are equal.
    def __lt__(self, other):
        if self.feature_id < other.feature_id:
            return True
        elif self.feature_id == other.feature_id:
            if self.key < other.key:
                return True
            elif self.key == other.key:
                if self.value < other.value:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

#this function takes a list of annotations and prints it as a 3-col table to a file
def write_annotations(annotations, file_out):
    # First, resolve duplicate gene names
    dups = {}
    for annotation in annotations:
        if annotation.key != "name":
            continue
        if annotation.value in dups:
            dups[annotation.value].append(annotation)
        else:
            dups[annotation.value] = [annotation]
    for dup in dups.values():
        if len(dup) > 1:
            for i, annotation in enumerate(dup):
                annotation.value += "_"+str(i)
    for annotation in annotations:
        file_out.write(annotation.feature_id+"\t"+annotation.key+"\t"+annotation.value+"\n") #write each annotation in the form: "feature_id, key, value"
