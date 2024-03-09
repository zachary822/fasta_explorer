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
import { useState } from "react";

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
        if (resp.status === 400) {
          form.setError("file", {
            type: "Bad Request",
            message: "Invalid File Type",
          });

          throw Error("Bad request");
        }

        return resp.json();
      })
      .then((data) => {
        setResult(data);
      });
  };

  return (
    <div>
      <Form {...form}>
        <form onSubmit={form.handleSubmit(submit)}>
          <FormField
            control={form.control}
            name="file"
            rules={{ required: true }}
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
          <Button type="submit" className="mt-8">
            Upload
          </Button>
        </form>
      </Form>
      <div className="text-base max-w-prose">
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
    </div>
  );
}
