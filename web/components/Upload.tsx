"use client";

import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import {
  Form,
  FormControl,
  FormDescription,
  FormField,
  FormItem,
  FormLabel,
  FormMessage,
} from "./ui/form";
import { FieldValues, SubmitHandler, useForm } from "react-hook-form";
import Plot from "react-plotly.js";
import { useState } from "react";
import { Textarea } from "./ui/textarea";

type FileResult = {
  sequences: {
    sequence: string;
    reverse_complement: string;
    gc_fraction: number;
  }[];
};

export default function Upload() {
  const form = useForm();

  const [result, setResult] = useState<FileResult>();

  const submit: SubmitHandler<FieldValues> = (_data, event) => {
    const body = new FormData(event?.target);

    fetch("/api/upload", {
      method: "POST",
      body,
    })
      .then((resp) => {
        if (resp.status >= 400) {
          resp.json().then((data) => {
            if (data.detail) {
              form.setError("file", {
                type: "Bad Request",
                message: data.detail,
              });
            } else if (
              Array.isArray(data) &&
              data.length > 0 &&
              data[0].type &&
              data[0].msg
            ) {
              form.setError("file", {
                type: data[0].type,
                message: data[0].msg,
              });
            }
          });
          return;
        }

        return resp.json();
      })
      .then((data) => {
        setResult(data);
      });
  };

  return (
    <>
      <div className="mb-4 w-full">
        <Form {...form}>
          <form
            onSubmit={form.handleSubmit(submit)}
            onReset={() => form.reset()}
          >
            <FormField
              control={form.control}
              name="file"
              render={({ field }) => (
                <FormItem>
                  <FormLabel>File</FormLabel>
                  <FormControl>
                    <Input type="file" accept=".fasta,.bam" {...field} />
                  </FormControl>
                  <FormDescription>Upload a FASTA or BAM file.</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />
            <FormField
              control={form.control}
              name="fasta"
              render={({ field }) => (
                <FormItem>
                  <FormLabel>FASTA</FormLabel>
                  <FormControl>
                    <Textarea {...field} />
                  </FormControl>
                  <FormDescription>
                    Input a fast formatted sequence.
                  </FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />
            <div className="mt-8">
              <Button type="submit" className="mr-4">
                Upload
              </Button>
              <Button type="reset" variant="destructive">
                Clear
              </Button>
            </div>
          </form>
        </Form>
      </div>
      {result && result.sequences.length > 1 && (
        <div>
          <Plot
            data={[
              {
                x: result.sequences.map((seq) => seq.sequence.length),
                type: "histogram",
              },
            ]}
            layout={{
              yaxis: { title: "count" },
              xaxis: { title: "length (bp)" },
            }}
          />
          <div className="text-base font-medium">
            Reads: <span>{result.sequences.length}</span>
          </div>
        </div>
      )}
      <div className="text-base w-full">
        {result &&
          result.sequences.map((seq, i) => (
            <div className="my-4" key={i}>
              <div className="text-base font-medium">
                Sequence <span>(length: {seq.sequence.length}bp)</span>
              </div>
              <div className="bg-white rounded-lg font-xs overflow-scroll p-4">
                {seq.sequence}
              </div>
              <div className="text-base font-medium">Reverse Complement</div>
              <div className="bg-white rounded-lg font-xs overflow-scroll p-4">
                {seq.reverse_complement}
              </div>
              <div className="text-base font-medium">
                GC Content:{" "}
                <span className="text-base font-mono">
                  {seq.gc_fraction.toLocaleString(undefined, {
                    style: "percent",
                    minimumFractionDigits: 2,
                  })}
                </span>
              </div>
            </div>
          ))}
      </div>
    </>
  );
}
